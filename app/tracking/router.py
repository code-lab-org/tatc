from celery import chain, group
from celery.result import GroupResult
from datetime import timezone
from fastapi import APIRouter, HTTPException
from geojson_pydantic import FeatureCollection
import json
import numpy as np
import pandas as pd
from uuid import UUID
from tatc.schemas.satellite import Satellite

from .schemas import GroundTrackRequest, GroundTrackSwathRequest
from .tasks import collect_ground_track_task, collect_ground_track_swath_task
from ..generation.schemas import TimeGenerator
from ..utils.schemas import CeleryTask, CeleryResultStatus, CeleryResultProgress
from ..utils.tasks import merge_feature_collections_task
from ..worker import app as celery_app

router = APIRouter()


@router.post("/orbit", response_model=CeleryTask)
async def collect_ground_track(request: GroundTrackRequest):
    if isinstance(request.times, TimeGenerator):
        times = pd.date_range(
            request.times.start, request.times.end, freq=request.times.delta
        ).tz_convert(timezone.utc)
    else:
        times = pd.DatetimeIndex(request.times).tz_convert(timezone.utc)
    task = chain(
        group(
            collect_ground_track_task.s(
                satellite.json(),
                satellite.instruments[request.instrument].json()
                if isinstance(request.instrument, int)
                and len(satellite.instruments) > request.instrument
                else request.instrument.json()
                if request.instrument is not None
                else None,
                [time.isoformat() for time in times],
                request.dict().get("mask", None),
            )
            for times in np.array_split(times, max(1, len(times) // 100))
            for satellite in request.satellite.generate_members()
        ),
        merge_feature_collections_task.s(),
    )()
    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {"task_id": task.id, "group_id": task.parent.id}
    else:
        return {"task_id": task.id}


@router.get("/orbit/{task_id}/status", response_model=CeleryResultStatus)
async def get_collect_ground_track_task_status(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"ready": task.ready(), "successful": task.successful()}


@router.get("/orbit/{group_id}/progress", response_model=CeleryResultProgress)
async def get_collect_ground_track_task_progress(group_id: UUID):
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"task_count": len(task.children), "completed_count": task.completed_count()}


@router.get(
    "/orbit/{task_id}/results",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={"200": {"content": {"application/x-msgpack": {}}}},
)
async def get_collect_ground_track_task_results(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    results = task.get()
    return FeatureCollection(features=json.loads(results).get("features"))


@router.post("/footprint", response_model=CeleryTask)
async def collect_ground_track_swath(request: GroundTrackSwathRequest):
    if isinstance(request.times, TimeGenerator):
        times = pd.date_range(
            request.times.start, request.times.end, freq=request.times.delta
        ).tz_convert(timezone.utc)
    else:
        times = pd.DatetimeIndex(request.times).tz_convert(timezone.utc)
    task = chain(
        group(
            collect_ground_track_swath_task.s(
                satellite.json(),
                satellite.instruments[request.instrument].json()
                if isinstance(request.instrument, int)
                and len(satellite.instruments) > request.instrument
                else request.instrument.json()
                if request.instrument is not None
                else None,
                [time.isoformat() for time in times],
                request.dict().get("mask", None),
                request.fast,
                request.resolution,
            )
            for times in np.array_split(
                times, max(1, len(times) // (100 if request.fast else 20))
            )
            for satellite in request.satellite.generate_members()
        ),
        body=merge_feature_collections_task.s(),
    )()
    if isinstance(task.parent, GroupResult):
        task.parent.save()
        return {"task_id": task.id, "group_id": task.parent.id}
    else:
        return {"task_id": task.id}


@router.get("/footprint/{task_id}/status", response_model=CeleryResultStatus)
async def get_collect_ground_track_swath_task_status(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {
        "ready": task.ready(),
        "successful": task.successful(),
        "completed_count": task.completed_count(),
    }


@router.get("/footprint/{group_id}/progress", response_model=CeleryResultProgress)
async def get_collect_ground_track_swath_task_progress(group_id: UUID):
    task = GroupResult.restore(str(group_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    return {"task_count": len(task.children), "completed_count": task.completed_count()}


@router.get(
    "/footprint/{task_id}/results",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={"200": {"content": {"application/x-msgpack": {}}}},
)
async def get_collect_ground_track_swath_task_results(task_id: UUID):
    task = celery_app.AsyncResult(str(task_id))
    if task is None:
        raise HTTPException(status_code=404, detail="Task not found.")
    if not task.ready():
        raise HTTPException(status_code=409, detail="Results not ready.")
    results = task.get()
    return FeatureCollection(features=json.loads(results).get("features"))
