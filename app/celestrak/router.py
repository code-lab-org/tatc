from fastapi import APIRouter
from fastapi.responses import PlainTextResponse
import aiohttp

router = APIRouter()


@router.get("/tle", response_class=PlainTextResponse)
async def get_celestrak_tle(name: str):
    async with aiohttp.ClientSession() as session:
        async with session.get(
            f"https://celestrak.com/NORAD/elements/gp.php?NAME={name}&FORMAT=TLE"
        ) as response:
            return await response.text()
