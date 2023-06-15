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
