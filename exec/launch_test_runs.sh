Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 2 &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f3 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f3 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier --fourier-levels 3 &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201030 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-stwa &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201030 --holdout 0 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-stwa &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 1 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 2 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 3 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 4 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex male --data-version 20201026 --model-version 20201029f2 --holdout 5 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 1 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 2 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 3 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 4 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "\n\n\n\n\n" &&
Rscript --vanilla ~/repos/covidemr/exec/03_ita_estimation_run.R --run-sex female --data-version 20201026 --model-version 20201029f2 --holdout 5 --use-covs intercept tfr unemp socserv tax_brackets hc_access elevation temperature --use-Z-sta --use-Z-fourier &&
printf "=============================================\n" &&
printf "==================== FIN ====================\n" &&
printf "=============================================\n";
