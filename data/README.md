# Data

## What is in this repository

`data/processed` holds the derived results used in the chapter:

- `pastas_trend_drought_summary.csv`: model fit, response times, trends, SGI drought statistics and latest status per monitoring station
- `drought_thresholds_monthly.csv`: monthly P50/P80/P90 depth thresholds per station
- `pastas_summary.csv`, `pastas_recharge_summary.csv`: earlier Pastas models
- `era5_pev_stations.csv`, `era5_tp_stations.csv`: daily ERA5-Land potential evaporation and precipitation at the stations (Copernicus)
- `borehole_national_stats.csv`, `district_borehole_summary.csv`, `corridor_district_summary.csv`, `region_summary.csv`: national, district and regional borehole statistics (medians and shares, no individual boreholes)
- `station_districts.csv`: district of each monitoring station and whether it lies in the cattle corridor

## What is not in this repository

The monitoring, borehole and drilling-log data belong to the **Ministry of Water and Environment (MWE), Uganda**, and are not redistributed here. They are available on request from the Directorate of Water Resources Management, MWE. To rerun the full workflow, place them as follows:

| Path | Content | Source |
|---|---|---|
| raw/Stations.csv | station register: name, number, coordinates, zone, purpose, data file | MWE |
| raw/levels/*.csv | daily depth to water (m bgl) | MWE monitoring network |
| raw/boreholes/BOREHOLE DATA.xlsx | national borehole database | MWE |
| raw/boreholes/Litho_strati.xlsx | Aswa drilling logs, sheet 'Aswa cat 1933-2012' | MWE |
| raw/era5/era5land_pev_t2m_YYYY.nc | ERA5-Land daily pev, tp, t2m, 1998–2026 | `scripts/get_era5_pev.py` (Copernicus CDS) |
| raw/dem/*.tif | 110 m Aswa and 220 m national DEM, EPSG:32636 | `scripts/get_dem.py` (Copernicus GLO-30) |
| raw/admin/uga_admin{0,1,2}.shp, Lakes.shp | admin boundaries, lakes | `scripts/get_admin_boundaries.py` (OCHA/UBOS via HDX) |

Borehole-level outputs (`boreholes_located.csv`, `aswa_layer_picks*.csv`) and the interactive 3D pages are produced by the scripts but are not published, because they list individual MWE boreholes.
