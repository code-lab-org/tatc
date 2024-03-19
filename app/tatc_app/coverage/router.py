# -*- coding: utf-8 -*-
"""
Router specifications for coverage analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from celery import group, chain
from celery.result import GroupResult
from fastapi import APIRouter, HTTPException
from geojson_pydantic import FeatureCollection
import json
from uuid import UUID

from .schemas import CoverageAnalysisRequest, CoverageAnalysisResult
from .tasks import run_coverage_analysis_task, grid_coverage_analysis_task
from ..generation.utils import generate_points, generate_cells
from ..celery.schemas import CeleryTask
from ..utils.tasks import merge_feature_collections_task
from ..worker import app as celery_app

router = APIRouter()


@router.post("/coverage", response_model=CeleryTask)
async def enqueue_coverage_analysis(request: CoverageAnalysisRequest):
    """
    Endpoint to enqueue a task to perform coverage analysis.

    Args:
        request (CoverageAnalysisRequest): User request.

    Returns:
        CeleryTask: Queued task.
    """
    task = chain(
        group(
            run_coverage_analysis_task.s(
                point.json(),
                [
                    satellite.json()
                    for constellation in request.satellites
                    for satellite in constellation.generate_members()
                ],
                request.start.isoformat(),
                request.end.isoformat(),
            )
            for point in generate_points(request.points)  # FIXME
        ),
        merge_feature_collections_task.s(),
        grid_coverage_analysis_task.s(
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
    "/coverage/{task_id}",
    response_model=CoverageAnalysisResult,
    response_model_exclude_none=True,
    responses={
        "200": {"content": {"application/json": {}, "application/x-msgpack": {}}}
    },
)
async def retrieve_coverage_analysis(task_id: UUID):
    """
    Endpoint to retrieve results of a completed coverage analysis task.

    Args:
        task_id (UUID): Task unique identifier.

    Returns:
        CoverageAnalysisResult: Coverage analysis results.
    """
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    return CoverageAnalysisResult.parse_raw(task.get())
