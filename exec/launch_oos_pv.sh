qsub -N f3fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3fal && sleep 0.25 &&
qsub -N f2fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2fal && sleep 0.25 &&
qsub -N f1fal_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1fal && sleep 0.25 &&
qsub -N f3fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3fal && sleep 0.25 &&
qsub -N f2fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2fal && sleep 0.25 &&
qsub -N f1fal_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1fal && sleep 0.25 &&
qsub -N f3fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3fal --oos && sleep 0.25 &&
qsub -N f2fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2fal --oos && sleep 0.25 &&
qsub -N f1fal_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1fal --oos && sleep 0.25 &&
qsub -N f3fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3fal --oos && sleep 0.25 &&
qsub -N f2fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2fal --oos && sleep 0.25 &&
qsub -N f1fal_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1fal --oos && sleep 0.25 &&
qsub -N f3fl_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3fl && sleep 0.25 &&
qsub -N f2fl_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2fl && sleep 0.25 &&
qsub -N f1fl_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1fl && sleep 0.25 &&
qsub -N f3fl_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3fl && sleep 0.25 &&
qsub -N f2fl_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2fl && sleep 0.25 &&
qsub -N f1fl_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1fl && sleep 0.25 &&
qsub -N f3fl_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3fl --oos && sleep 0.25 &&
qsub -N f2fl_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2fl --oos && sleep 0.25 &&
qsub -N f1fl_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1fl --oos && sleep 0.25 &&
qsub -N f3fl_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3fl --oos && sleep 0.25 &&
qsub -N f2fl_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2fl --oos && sleep 0.25 &&
qsub -N f1fl_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1fl --oos && sleep 0.25 &&
qsub -N f3_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3 && sleep 0.25 &&
qsub -N f2_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2 && sleep 0.25 &&
qsub -N f1_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1 && sleep 0.25 &&
qsub -N f3_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3 && sleep 0.25 &&
qsub -N f2_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2 && sleep 0.25 &&
qsub -N f1_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1 && sleep 0.25 &&
qsub -N f3_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f3 --oos && sleep 0.25 &&
qsub -N f2_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f2 --oos && sleep 0.25 &&
qsub -N f1_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113f1 --oos && sleep 0.25 &&
qsub -N f3_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f3 --oos && sleep 0.25 &&
qsub -N f2_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f2 --oos && sleep 0.25 &&
qsub -N f1_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113f1 --oos && sleep 0.25 &&
qsub -N nof_m_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113nof && sleep 0.25 &&
qsub -N nof_f_is ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113nof && sleep 0.25 &&
qsub -N nof_m_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210113nof --oos && sleep 0.25 &&
qsub -N nof_f_oos ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20210113 --model-version 20210113nof --oos && sleep 0.25 &&echo "================================================================" &&
echo "============================= DONE =============================" &&
echo "================================================================";
