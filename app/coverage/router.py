from celery import group, chain
from celery.result import GroupResult
from fastapi import APIRouter, HTTPException
from geojson_pydantic import FeatureCollection
import json
from uuid import UUID

from .schemas import CollectCoverageStatisticsRequest, CoverageStatisticsAnalysisResult
from .tasks import compute_coverage_statistics_task, aggregate_coverage_statistics_task
from ..generation.utils import generate_points, generate_cells
from ..utils.schemas import CeleryTask, CeleryResultStatus, CeleryResultProgress
from ..utils.tasks import merge_feature_collections_task
from ..worker import app as celery_app

router = APIRouter()


@router.post("/coverage", response_model=CeleryTask)
async def collect_coverage_statistics(request: CollectCoverageStatisticsRequest):
    task = chain(
        group(
            compute_coverage_statistics_task.s(
                point.json(),
                [
                    satellite.json()
                    for constellation in request.satellites
                    for satellite in constellation.generate_members()
                ],
                request.start.isoformat(),
                request.end.isoformat(),
                True,
                request.sample_distance,
            )
            for point in generate_points(request.points)
        ),
        merge_feature_collections_task.s(),
        aggregate_coverage_statistics_task.s(
            generate_cells(request.cells).to_json(show_bbox=False, drop_id=True)
        ),
    )()
    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {"task_id": task.id, "group_id": task.parent.id}
    else:
        return {"task_id": task.id}


@router.get("/coverage/{task_id}/status", response_model=CeleryResultStatus)
async def get_collect_coverage_statistics_task_status(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"ready": task.ready(), "successful": task.successful()}


@router.get("/coverage/{group_id}/progress", response_model=CeleryResultProgress)
async def get_collect_coverage_statistics_task_progress(group_id: UUID):
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"task_count": len(task.children), "completed_count": task.completed_count()}


@router.get(
    "/coverage/{task_id}/results",
    response_model=CoverageStatisticsAnalysisResult,
    response_model_exclude_none=True,
    responses={
        "200": {"content": {"application/json": {}, "application/x-msgpack": {}}}
    },
)
async def get_collect_coverage_statistics_task_results(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    return CoverageStatisticsAnalysisResult.parse_raw(task.get())
