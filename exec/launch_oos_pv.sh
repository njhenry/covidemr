qsub -N f3fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f3fal && sleep 0.25 &&
qsub -N f2fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f2fal && sleep 0.25 &&
qsub -N f1fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f1fal && sleep 0.25 &&
qsub -N nof_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213nof && sleep 0.25 &&
qsub -N f3fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f3fal && sleep 0.25 &&
qsub -N f2fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f2fal && sleep 0.25 &&
qsub -N f1fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f1fal && sleep 0.25 &&
qsub -N nof_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213nof && sleep 0.25 &&
qsub -N f3fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f3fal --oos && sleep 0.25 &&
qsub -N f2fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f2fal --oos && sleep 0.25 &&
qsub -N f1fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213f1fal --oos && sleep 0.25 &&
qsub -N nof_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201213 --model-version 20201213nof --oos && sleep 0.25 &&
qsub -N f3fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f3fal --oos && sleep 0.25 &&
qsub -N f2fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f2fal --oos && sleep 0.25 &&
qsub -N f1fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213f1fal --oos && sleep 0.25 &&
qsub -N nof_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201213 --model-version 20201213nof --oos && sleep 0.25 &&echo "================================================================" &&
echo "============================= DONE =============================" &&
echo "================================================================";
