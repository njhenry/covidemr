#!/bin/bash
#! -q short.qc
#! -pe shmem 10
#! -e /users/moore/bsb907/logging/
#! -o /users/moore/bsb907/logging/
#! -m n
#! -M nathaniel.henry@ndm.ox.ac.uk
# R submission script
/users/moore/bsb907/miniconda3/envs/r_env/bin/Rscript --vanilla "$@"
