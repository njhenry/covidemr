qsub -N is_m_f3 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f3fageloc && sleep 1 &&
qsub -N is_m_f2 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f2fageloc && sleep 1 &&
qsub -N is_m_f1 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f1fageloc && sleep 1 &&
qsub -N is_f_f3 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f3fageloc && sleep 1 &&
qsub -N is_f_f2 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f2fageloc && sleep 1 &&
qsub -N is_f_f1 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f1fageloc && sleep 1 &&
qsub -N oos_m_f3 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f3fageloc --oos && sleep 1 &&
qsub -N oos_m_f2 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f2fageloc --oos && sleep 1 &&
qsub -N oos_m_f1 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex male --data-version 20201203 --model-version 20201211f1fageloc --oos && sleep 1 &&
qsub -N oos_f_f3 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f3fageloc --oos && sleep 1 &&
qsub -N oos_f_f2 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f2fageloc --oos && sleep 1 &&
qsub -N oos_f_f1 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/04b_predictive_validity_metrics.R --run-sex female --data-version 20201203 --model-version 20201211f1fageloc --oos && sleep 1 &&
echo "================================================================" &&
echo "============================= DONE =============================" &&
echo "================================================================";
