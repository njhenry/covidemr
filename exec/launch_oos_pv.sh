qsub -N f3fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f3fal && sleep 0.25 &&
qsub -N f2fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f2fal && sleep 0.25 &&
qsub -N f1fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f1fal && sleep 0.25 &&
qsub -N nof_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212nof && sleep 0.25 &&
qsub -N f3fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f3fal && sleep 0.25 &&
qsub -N f2fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f2fal && sleep 0.25 &&
qsub -N f1fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f1fal && sleep 0.25 &&
qsub -N nof_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212nof && sleep 0.25 &&
qsub -N f3fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f3fal --oos && sleep 0.25 &&
qsub -N f2fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f2fal --oos && sleep 0.25 &&
qsub -N f1fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212f1fal --oos && sleep 0.25 &&
qsub -N nof_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201212 --model-version 20201212nof --oos && sleep 0.25 &&
qsub -N f3fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f3fal --oos && sleep 0.25 &&
qsub -N f2fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f2fal --oos && sleep 0.25 &&
qsub -N f1fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212f1fal --oos && sleep 0.25 &&
qsub -N nof_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201212 --model-version 20201212nof --oos && sleep 0.25 &&echo "================================================================" &&
echo "============================= DONE =============================" &&
echo "================================================================";
