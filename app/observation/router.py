from celery import chain, group
from celery.result import GroupResult
from fastapi import APIRouter, HTTPException
from geojson_pydantic import FeatureCollection
import json
from uuid import UUID

from .schemas import CollectObservationsRequest
from .tasks import collect_observations_task, aggregate_observations_task
from ..generation.utils import generate_points, generate_cells
from ..utils.schemas import CeleryTask, CeleryResultStatus, CeleryResultProgress
from ..utils.tasks import merge_feature_collections_task
from ..worker import app as celery_app

router = APIRouter()


@router.post(
    "/observations", response_model=CeleryTask, response_model_exclude_none=True
)
async def create_collect_observations_task(request: CollectObservationsRequest):
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
                sample_distance=request.sample_distance
            )
            for point in generate_points(request.points)
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


@router.get("/observations/{task_id}/status", response_model=CeleryResultStatus)
async def get_collect_observations_task_status(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"ready": task.ready(), "successful": task.successful()}


@router.get("/observations/{group_id}/progress", response_model=CeleryResultProgress)
async def get_collect_observations_task_progress(group_id: UUID):
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"task_count": len(task.children), "completed_count": task.completed_count()}


@router.get(
    "/observations/{task_id}/results",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={
        "200": {"content": {"application/json": {}, "application/x-msgpack": {}}}
    },
)
async def get_collect_observations_task_results(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    results = task.get()
    return FeatureCollection(features=json.loads(results).get("features"))
