qsub -N f_f2fal_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f2fal_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N m_f2fal_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f2fal_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N f_f2fal_ns_nn ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f2fal_ns_nn --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --fourier-ns && sleep 0.5 &&
qsub -N m_f2fal_ns_nn ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f2fal_ns_nn --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --fourier-ns && sleep 0.5 &&
qsub -N f_f3fal_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f3fal_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 --fourier-groups age_group_code location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N m_f3fal_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f3fal_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 --fourier-groups age_group_code location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N f_f2fl_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f2fl_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N m_f2fl_ns ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f2fl_ns --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups location_code --use-nugget --fourier-ns && sleep 0.5 &&
qsub -N f_f2fal_s ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f2fal --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --use-nugget && sleep 0.5 &&
qsub -N m_f2fal_s ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f2fal --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code --use-nugget && sleep 0.5 &&
qsub -N f_f2fal_s_nn ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20210113 --model-version 20210421_f2fal_nn --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code && sleep 0.5 &&
qsub -N m_f2fal_s_nn ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20210113 --model-version 20210421_f2fal_nn --holdout 0 --use-covs intercept year_cov tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 --fourier-groups age_group_code location_code && sleep 0.5;
