conda activate base &&
  conda create -n r_env -y -c conda-forge r-base=4.0.3 r-essentials r-renv r-devtools &&
  conda activate r_env;
