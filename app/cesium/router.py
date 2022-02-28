from fastapi import APIRouter
from fastapi.responses import PlainTextResponse

def get_cesium_router(token):
    router = APIRouter()

    @router.get("/token", response_class=PlainTextResponse)
    async def get_cesium_token():
        return token

    return router
