# ccmmma-postpro

Post-processing utilities for coastal and atmospheric model outputs produced at [Meteo@Uniparthenope](https://meteo.uniparthenope.it). The scripts read NetCDF outputs from WRF, ROMS, WaveWatch III, and WACOMM, interpolate them onto user-friendly grids, and export compact archive files for downstream visualization or delivery.

## Requirements
- Python 3.11 (refer to [wrf-python](https://wrf-python.readthedocs.io) for compatibility limitations.
- NetCDF4-compatible build environment
- Dependencies listed in `requirements.txt`

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Usage
Each script expects specific inputs:

- **WRF:** `python postpro-wrf5.py <init_date> <current.nc> <previous_hour.nc> <midnight.nc> <output.nc>`
- **ROMS:** `python postpro-rms3.py <init_date> <source.nc> <history_dir> <output.nc>`
- **WaveWatch III:** `python postpro-ww33.py <init_date> <source.nc> <history_dir> <output.nc>`
- **WACOMM:** `python postpro-wcm3.py <init_date> <source.nc> <history_dir> <output.nc>`

Logs are emitted at `INFO` level to describe each processing step, making it easier to follow how fields are interpolated or derived.

## Notes
- The code assumes inputs follow the grid conventions used by Meteo@Uniparthenope operational systems; adjust file paths to match your local archive layout.
- Historical files (previous hour and midnight runs) enable the scripts to compute deltas such as hourly rain or wind shifts.
