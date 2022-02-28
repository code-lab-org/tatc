from pydantic import BaseModel, Field
from typing import Optional, List
from uuid import UUID


class CeleryTask(BaseModel):
    task_id: UUID = Field(..., description="Unique task identifier.")
    group_id: Optional[UUID] = Field(
        None, description="Unique group task identifier (if applicable)."
    )


class CeleryResultProgress(BaseModel):
    task_count: int = Field(..., description="Number of sub-tasks to complete.")
    completed_count: int = Field(..., description="Number of completed sub-tasks.")


class CeleryResultStatus(BaseModel):
    ready: bool = Field(..., description="True, if this task is ready.")
    successful: bool = Field(..., description="True, if this task was successful.")
