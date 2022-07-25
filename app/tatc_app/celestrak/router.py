from fastapi import APIRouter
from fastapi.responses import PlainTextResponse
import aiohttp


router = APIRouter()


@router.get("/tle", response_class=PlainTextResponse)
async def get_celestrak_tle(name: str):
    """
    Endpoint to make general perturbations (GP) requests to CelesTrak.
    Allows clients to make requests without cross-site scripting requests.

    Args:
        name (str): the satellite name

    Returns:
        PlainTextResponse: the GP in two line elements (TLE) format
    """
    async with aiohttp.ClientSession() as session:
        async with session.get(
            f"https://celestrak.com/NORAD/elements/gp.php?NAME={name}&FORMAT=TLE"
        ) as response:
            result = await response.text()
    return result
