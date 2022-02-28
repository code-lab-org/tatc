from celery import (
    chain,
    group,
    chord
)
from celery.result import GroupResult
from fastapi import (
    APIRouter,
    HTTPException
)
from geojson_pydantic import FeatureCollection
import json
from uuid import UUID

from ..worker import app as celery_app
from .tasks import (
    collect_downlinks_task,
    compute_latency_task,
    aggregate_downlinks_task,
    aggregate_latencies_task
)
from ..utils.schemas import (
    CeleryTask,
    CeleryResultStatus,
    CeleryResultProgress
)
from .schemas import ComputeLatenciesRequest
from ..generation.utils import generate_points

router = APIRouter()


@router.post('/latencies',
    response_model=CeleryTask
)
async def create_compute_latencies_task(request: ComputeLatenciesRequest):

    task = chain(
                chord(
                    group([
                        collect_downlinks_task.s(
                            station.json(),
                            satellite.json(),
                            request.start.isoformat(),
                            request.end.isoformat()
                        )
                        for station in request.stations for satellite in request.satellites
                    ]),
                    body=aggregate_downlinks_task.s()
                    ),
                chord(
                    group([
                        compute_latency_task.s(
                            point.json(),
                            satellite.json(),
                            instrument.json(),
                            request.start.isoformat(),
                            request.end.isoformat(),
                            request.omit_solar
                            )
                        for point in generate_points(request.points) for satellite in request.satellites for instrument in satellite.instruments
                    ]),
                    body=aggregate_latencies_task.s()
                )
            )()


    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {
            "task_id": task.id,
            "group_id": task.parent.id
        }
    else:
        return {
            "task_id": task.id
        }


@router.get('/latencies/{task_id}/status', response_model=CeleryResultStatus)
async def get_compute_latencies_task_status(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"ready": task.ready(), "successful": task.successful()}

@router.get(
    '/latencies/{group_id}/progress',
    response_model=CeleryResultProgress
)
async def get_collect_latencies_task_progress(group_id: UUID):
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {
        "task_count": len(task.children),
        "completed_count": task.completed_count()
    }


@router.get(
    '/latencies/{task_id}/results',
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={
        '200': {
            'content': {
                'application/json': { },
                'application/x-msgpack': { }
            }
        }
    }
)
async def get_compute_latencies_task_results(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    results = task.get()
    return FeatureCollection(
        features=json.loads(results).get('features')
    )
