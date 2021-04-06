##
## Execute this script with a single argument to specify the environment name
## Default environment is "r_env"
## Example: ./conda_create_R_environment.sh my_new_env
##
ENV_NAME=${1:-r_env};
CRAN_REPO=https://cloud.r-project.org;
INLA_REPO=https://inla.r-inla-download.org/R/stable;

conda activate base &&
  conda create -n $ENV_NAME -y -c conda-forge r-base=4 gxx_linux-64 gdal \
    r-codetools r-sf r-devtools r-renv &&
  conda activate $ENV_NAME &&
  pip install -U radian &&
  R --vanilla -e "devtools::install_github('jalvesaq/colorout')" &&
  R --vanilla -e "install.packages('INLA',repos=c(CRAN='$CRAN_REPO',INLA='$INLA_REPO'),dep=TRUE)";
