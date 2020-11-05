Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201104 --model-version 20201104f3fageloc --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 --fourier-groups age_group_code location_code && printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201104 --model-version 20201104f3fageloc --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 --fourier-groups age_group_code location_code && printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201104 --model-version 20201104f3fageloc &&
Rscript --vanilla ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201104 --model-version 20201104f3fageloc &&
echo "================================================================" &&
echo "============================= DONE =============================" &&
echo "================================================================";
