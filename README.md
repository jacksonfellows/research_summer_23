# Python environment  #

To build the lock file:
```
conda-lock -f environment.yml
```
To create a Conda environment from the lock file:
```
conda-lock install --name YOUR_ENVIRONMENT_NAME_HERE conda-lock.yml
```

# Stuff you need to download #

Can download USGS quaternary faults here: <https://earthquake.usgs.gov/static/lfs/nshm/qfaults/Qfaults_GIS.zip>.
