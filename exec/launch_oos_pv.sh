qsub -N m_f2fal ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20210113 --model-version 20210421_f2fal && sleep 0.25 &&
echo "============================= DONE =============================" &&
echo "================================================================";
