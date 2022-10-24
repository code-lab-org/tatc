from fastapi import APIRouter, HTTPException
from uuid import UUID
from celery.result import GroupResult

from .schemas import CeleryResultStatus, CeleryResultProgress
from ..worker import app as celery_app


router = APIRouter()


@router.get("/{task_id}/status", response_model=CeleryResultStatus)
async def get_collect_coverage_statistics_task_status(task_id: UUID):
    """
    Endpoint to retrieve enqueued task status.

    Args:
        task_id (UUID): Task unique identifier.

    Returns:
        CeleryResultStatus: Task results status.
    """
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"ready": task.ready(), "successful": task.successful()}


@router.get("/{group_id}/progress", response_model=CeleryResultProgress)
async def get_collect_coverage_statistics_task_progress(group_id: UUID):
    """
    Endpoint to retrieve enqueued group task progress.

    Args:
        group_id (UUID): Group task unique identifier.

    Returns:
        CeleryResultProgress: Group task progress.
    """
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"task_count": len(task.children), "completed_count": task.completed_count()}
