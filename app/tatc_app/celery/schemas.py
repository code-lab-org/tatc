# -*- coding: utf-8 -*-
"""
Schema specifications for celery endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Optional
from uuid import UUID

from pydantic import BaseModel, Field


class CeleryTask(BaseModel):
    """
    Enqueued Celery task.
    """

    task_id: UUID = Field(..., description="Unique task identifier.")
    group_id: Optional[UUID] = Field(
        None, description="Unique group task identifier (if applicable)."
    )


class CeleryResultProgress(BaseModel):
    """
    Progress of an enqueued Celery group task.
    """

    task_count: int = Field(..., description="Number of sub-tasks to complete.")
    completed_count: int = Field(..., description="Number of completed sub-tasks.")


class CeleryResultStatus(BaseModel):
    """
    Status of an enqueued Celery task.
    """

    ready: bool = Field(..., description="True, if this task is ready.")
    successful: bool = Field(..., description="True, if this task was successful.")
