# Python environment  #

To build the lock file:
```
conda-lock -f environment.yml
```
To create a Conda environment from the lock file:
```
conda-lock install --name YOUR_ENVIRONMENT_NAME_HERE conda-lock.yml
```

# Stuff to download #

Can download USGS quaternary faults here: <https://earthquake.usgs.gov/static/lfs/nshm/qfaults/Qfaults_GIS.zip>.
Can find volcanoes for a search range here: <https://www.ngdc.noaa.gov/hazel/view/hazards/volcano/loc-search>.

# Testing #

To run all tests run:
```
pytest
```

For less output run:
```
pytest -q --disable-warnings
```

# Notes #

1964 rupture traced from Suleimani, E., Nicolsky, D. J., Haeussler, P. J., & Hansen, R. (2011). Combined Effects of Tectonic and Landslide-Generated Tsunami Runup at Seward, Alaska During the M W 9.2 1964 Earthquake. Pure and Applied Geophysics, 168(6–7), 1053–1074. https://doi.org/10.1007/s00024-010-0228-4
All other ruptures traced (and all rupture magnitudes taken) from Liu, C., Lay, T., & Xiong, X. (2022). The 29 July 2021 MW 8.2 Chignik, Alaska Peninsula Earthquake Rupture Inferred From Seismic and Geodetic Observations: Re-Rupture of the Western 2/3 of the 1938 Rupture Zone. Geophysical Research Letters, 49(4), e2021GL096004. https://doi.org/10.1029/2021GL096004
