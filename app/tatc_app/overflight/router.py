from celery import chain, group
from celery.result import GroupResult
from fastapi import APIRouter, HTTPException
from geojson_pydantic import FeatureCollection
import json
from uuid import UUID

from .schemas import OverflightAnalysisRequest
from .tasks import collect_observations_task, aggregate_observations_task
from ..generation.utils import generate_points, generate_cells
from ..celery.schemas import CeleryTask
from ..utils.tasks import merge_feature_collections_task
from ..worker import app as celery_app

router = APIRouter()


@router.post("/overflight", response_model=CeleryTask, response_model_exclude_none=True)
async def enqueue_overflight_analysis(request: OverflightAnalysisRequest):
    """
    Endpoint to enqueue a task to perform overflight analysis.

    Args:
        request (OverflightAnalysisRequest): User request.

    Returns:
        CeleryTask: Queued task.
    """
    task = chain(
        group(
            collect_observations_task.s(
                point.json(),
                satellite.json(),
                satellite.instruments[request.instrument].json()
                if isinstance(request.instrument, int)
                else request.instrument.json(),
                request.start.isoformat(),
                request.end.isoformat(),
                request.omit_solar,
            )
            for point in generate_points(request.points)  # FIXME
            for satellite in request.satellite.generate_members()
        ),
        merge_feature_collections_task.s(),
        aggregate_observations_task.s(),
    )()
    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {"task_id": task.id, "group_id": task.parent.id}
    else:
        return {"task_id": task.id}


@router.get(
    "/overflight/{task_id}",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={
        "200": {"content": {"application/json": {}, "application/x-msgpack": {}}}
    },
)
async def retrieve_overflight_anlaysis(task_id: UUID):
    """
    Endpoint to retrieve overflight analysis results.

    Args:
        task_id (UUID): Task unique identifier.

    Returns:
        FeatureCollection: Overflight analysis results.
    """
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    results = task.get()
    return FeatureCollection(features=json.loads(results).get("features"))
