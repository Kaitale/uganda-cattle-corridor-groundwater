@echo off
title Groundwater chapter - full workflow (leave this window open)
cd /d "%~dp0"
python -c "import pandas" >nul 2>&1
if errorlevel 1 (
  for %%C in ("%USERPROFILE%\anaconda3\Scripts\activate.bat" "%USERPROFILE%\miniconda3\Scripts\activate.bat" "C:\ProgramData\anaconda3\Scripts\activate.bat") do if exist %%C call %%C >nul 2>&1
)
echo Checking packages ...
python -m pip install --quiet -r requirements.txt
if not exist data\raw\admin\uga_admin2.shp python scripts\get_admin_boundaries.py "F:\APPs\Semuliki_HEC_RTS\Data\uga_admin_boundaries.shp"
if not exist data\raw\admin\Lakes.shp copy "F:\APPs\KaiEWS\Data\River gauge data\WL_gauge_LocationMap\Lakes.*" data\raw\admin\ >nul
echo.
echo Running workflow (add --era5 in this file to re-download ERA5-Land) ...
python scripts\run_all.py > logs_run_all.txt 2>&1
type logs_run_all.txt
echo.
echo Finished. Tables: data\processed   Figures: figures\paper
pause
