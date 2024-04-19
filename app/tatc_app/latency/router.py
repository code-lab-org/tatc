# -*- coding: utf-8 -*-
"""
Router specifications for latency analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from uuid import UUID

from celery import chain, group
from celery.result import GroupResult
from fastapi import APIRouter, HTTPException

from ..worker import app as celery_app
from .tasks import (
    collect_downlinks_task,
    run_latency_analysis_task,
    grid_latency_analysis_task,
)
from ..celery.schemas import CeleryTask
from ..utils.tasks import merge_feature_collections_task
from .schemas import LatencyAnalysisRequest, LatencyAnalysisResult
from ..generation.utils import generate_points, generate_cells

router = APIRouter()


@router.post("/latency", response_model=CeleryTask)
async def enqueue_latency_analysis(request: LatencyAnalysisRequest):
    """
    Endpoint to enqueue a task to perform latency analysis.

    Args:
        request (LatencyAnalysisRequest): User request.

    Returns:
        CeleryTask: Queued task.
    """
    task = chain(
        group(
            collect_downlinks_task.s(
                [station.model_dump_json() for station in request.stations],
                satellite.model_dump_json(),
                request.start.isoformat(),
                request.end.isoformat(),
            )
            for constellation in request.satellites
            for satellite in constellation.generate_members()
        ),
        merge_feature_collections_task.s(),
        group(
            run_latency_analysis_task.s(
                point.model_dump_json(),
                [
                    satellite.model_dump_json()
                    for constellation in request.satellites
                    for satellite in constellation.generate_members()
                ],
                request.start.isoformat(),
                request.end.isoformat(),
            )
            for point in generate_points(request.points)  # FIXME
        ),
        merge_feature_collections_task.s(),
        grid_latency_analysis_task.s(
            generate_cells(request.cells).to_json(
                show_bbox=False, drop_id=True
            )  # FIXME
        ),
    )()
    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {"task_id": task.id, "group_id": task.parent.id}
    else:
        return {"task_id": task.id}


@router.get(
    "/latency/{task_id}",
    response_model=LatencyAnalysisResult,
    response_model_exclude_none=True,
    responses={
        "200": {"content": {"application/json": {}, "application/x-msgpack": {}}}
    },
)
async def retrieve_latency_analysis(task_id: UUID):
    """
    Endpoint to retrieve results of a completed latency analysis task.

    Args:
        task_id (UUID): Task unique identifier.

    Returns:
        LatencyAnalysisResult: Latency analysis results.
    """
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    return LatencyAnalysisResult.model_validate_json(task.get())
