qsub viz_f3fal ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f3fal && sleep 0.25 &&
qsub viz_f2fal ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f2fal && sleep 0.25 &&
qsub viz_f1fal ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f1fal && sleep 0.25 &&
qsub viz_f3fl ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f3fl && sleep 0.25 &&
qsub viz_f2fl ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f2fl && sleep 0.25 &&
qsub viz_f1fl ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f1fl && sleep 0.25 &&
qsub viz_f3 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f3 && sleep 0.25 &&
qsub viz_f2 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f2 && sleep 0.25 &&
qsub viz_f1 ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214f1 && sleep 0.25 &&
qsub viz_nof ~/repos/covidemr/exec/cluster_R_exec.sh ~/repos/covidemr/exec/viz/map_excess_mort_weekly.R --data-version 20201214 --model-version 20201214nof && sleep 0.25;