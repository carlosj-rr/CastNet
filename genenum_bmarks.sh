#mkdir 1-30_GeneNum_Bmarks
cd ~/Programs/CastNet/1-30_GeneNum_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
for rep in {1..30}
do
  cp CastNet_parameters.py CastNet_parameters_${rep}.py
  cp CastNet.py CastNet_${rep}.py
  old_rate=$(grep seq_mutation_rate CastNet_parameters_${rep}.py | tr [:blank:] "," | cut -d"," -f3)
  old_num_genes=$(grep num_genes CastNet_parameters_${rep}.py | tr [:blank:] "," | cut -d"," -f3)
  new_rate=$(sed -n ${rep}p evolrates.csv)
  sed -i "s/$old_rate/$new_rate/" CastNet_parameters_${rep}.py
  sed -i "s/num_genes = $old_num_genes/num_genes = ${rep}/" CastNet_parameters_${rep}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${rep}/" CastNet_${rep}.py
  beginning=$(date +%H:%M:%S)
  python CastNet_${rep}.py ${rep}_genesRun
  end=$(date +%H:%M:%S)
  echo $rep,$beginning,$end >> muts_parallel_bmark_2.csv
done

echo "Finished first bit successfully."
