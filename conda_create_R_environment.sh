##
## Execute this script with a single argument to specify the environment name
## Default environment is "r_env"
## Example: ./conda_create_R_environment.sh my_new_env
##
ENV_NAME=${1:-r_env}

conda activate base &&
  conda create -n $ENV_NAME -y -c conda-forge r-base=4 gxx_linux-64 gdal r-sf \
    r-devtools r-renv &&
  conda activate $ENV_NAME &&
  pip install -U radian &&
  R --vanilla -e "devtools::install_github('jalvesaq/colorout')";
