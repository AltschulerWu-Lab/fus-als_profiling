declare -a MODEL=("irf" "fcnn" "logistic")
declare -a CYTOSOLIC=("TRUE" "FALSE")

for m in "${MODEL[@]}"
do
  
  for c in "${CYTOSOLIC[@]}" 
  do
    echo "RUNNING: $m, $c"
    Rscript train_fus_models.R --MODEL "$m" --CYTOSOLIC "$c"
    Rscript fig_scores.R --MODEL "$m" --CYTOSOLIC "$c"
    echo "####################"
  done

done

