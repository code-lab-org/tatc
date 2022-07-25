from fastapi import APIRouter
from fastapi.responses import PlainTextResponse


def get_cesium_router(token: str) -> APIRouter:
    """
    Router that provides a single endpoint to retrieve the Cesium access token.

    Args:
        token (str): the Cesium access token

    Returns:
        APIRouter: the FastAPI router
    """
    router = APIRouter()

    @router.get("/token", response_class=PlainTextResponse)
    async def get_cesium_token():
        """
        Endpoint to retrieve the Cesium access token.

        Returns:
            PlainTextResponse: the response containing the access token
        """
        return token

    return router
