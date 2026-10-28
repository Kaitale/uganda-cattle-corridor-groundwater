"""Download Copernicus GLO-30 DEM for the Aswa catchment / northern Uganda from OpenTopography.
Run on a machine with internet. Free API key: https://portal.opentopography.org (My Account > API key)."""
import requests, sys
KEY = sys.argv[1] if len(sys.argv) > 1 else "PASTE_YOUR_KEY"
url = "https://portal.opentopography.org/API/globaldem"
p = dict(demtype="COP30", south=1.5, north=4.3, west=31.3, east=34.6, outputFormat="GTiff", API_Key=KEY)
r = requests.get(url, params=p, timeout=600); r.raise_for_status()
open("../data/raw/dem/cop30_north_uganda.tif", "wb").write(r.content)
print("saved", len(r.content) / 1e6, "MB")
