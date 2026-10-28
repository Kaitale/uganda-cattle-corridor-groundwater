"""Get Uganda administrative boundaries (OCHA/UBOS COD-AB, levels 0-2) into data/raw/admin.

Uses a local copy if you already have one (pass its folder), otherwise downloads the latest
shapefile zip from the Humanitarian Data Exchange (CC BY-IGO).
    python get_admin_boundaries.py                      # download from HDX
    python get_admin_boundaries.py "F:/APPs/Semuliki_HEC_RTS/Data/uga_admin_boundaries.shp"
A lakes layer (Lakes.shp, any CRS) should also be placed in data/raw/admin for the maps.
"""
import glob, io, os, shutil, sys, zipfile
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "..", "data", "raw", "admin")
HDX = "https://data.humdata.org/api/3/action/package_show?id=cod-ab-uga"


def from_folder(src):
    for f in glob.glob(os.path.join(src, "uga_admin[012].*")):
        shutil.copy(f, OUT); print("copied", os.path.basename(f))


def from_hdx():
    import requests
    res = requests.get(HDX, timeout=60).json()["result"]["resources"]
    url = next(r["url"] for r in res if r["url"].lower().endswith(".zip") and "shp" in r["url"].lower())
    print("downloading", url)
    z = zipfile.ZipFile(io.BytesIO(requests.get(url, timeout=300).content))
    for n in z.namelist():
        b = os.path.basename(n)
        if b.startswith(("uga_admin0", "uga_admin1", "uga_admin2")):
            open(os.path.join(OUT, b), "wb").write(z.read(n)); print("saved", b)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    from_folder(sys.argv[1]) if len(sys.argv) > 1 else from_hdx()
