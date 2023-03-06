cd ~/Programs/CastNet/1-20_GeneNum_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
for rep in {1..20}
do
  cp CastNet_parameters.py CastNet_parameters_${rep}.py
  cp CastNet.py CastNet_${rep}.py
  old_rate=$(grep seq_mutation_rate CastNet_parameters_${rep}.py | tr [:blank:] "," | cut -d"," -f3)
  old_num_genes=$(grep num_genes CastNet_parameters_${rep}.py | tr [:blank:] "," | cut -d"," -f3)
  new_rate=$(sed -n ${rep}p evolrates.csv)
  sed -i "s/$old_rate/$new_rate/" CastNet_parameters_${rep}.py
  sed -i "s/num_genes = $old_num_genes/num_genes = ${rep}/" CastNet_parameters_${rep}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${rep}/" CastNet_${rep}.py
  date +%H,%M,%S > times_${rep}-genes.csv
  python CastNet_${rep}.py ${rep}_genesRun
  date +%H,%M,%S >> times_${rep}-genes.csv
done

echo "Finished first bit successfully."

cd ~/Programs/CastNet/10-1K_PopSize_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
for line in {1..8}
do
  new_pop_size=$(sed -n ${line}p sizes_n_rates.csv | cut -d"," -f1)
  new_rate=$(sed -n ${line}p sizes_n_rates.csv | cut -d"," -f2)
  cp CastNet_parameters.py CastNet_parameters_${new_pop_size}.py
  cp CastNet.py CastNet_${new_pop_size}.py
  old_rate=$(grep seq_mutation_rate CastNet_parameters_${new_pop_size}.py | tr [:blank:] "," | cut -d"," -f3)
  old_pop_size=$(egrep -o "pop_size = [0-9]+" CastNet_parameters_${new_pop_size}.py)
  sed -i "s/$old_rate/$new_rate/" CastNet_parameters_${new_pop_size}.py
  sed -i "s/$old_pop_size/pop_size = $new_pop_size/" CastNet_parameters_${new_pop_size}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${new_pop_size}/" CastNet_${new_pop_size}.py
  date +%H,%M,%S > start_stop_${new_pop_size}-indivs.csv
  python CastNet_${new_pop_size}.py ${new_pop_size}_indivsRun
  date +%H,%M,%S >> start_stop_${new_pop_size}-indivs.csv
done
cd ~/Programs/CastNet
echo "Finished second bit successfully"

cd ~/Programs/CastNet/1e-8_1e-1_MutRate_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .

for line in {1..8}
do
  new_rate_name=$(sed -n ${line}p rates_table.list | cut -d"," -f1)
  new_rate=$(sed -n ${line}p rates_table.list | cut -d"," -f2)
  cp CastNet_parameters.py CastNet_parameters_${new_rate_name}.py
  cp CastNet.py CastNet_${new_rate_name}.py
  old_rate=old_rate=$(grep seq_mutation_rate CastNet_parameters_${new_rate_name}.py | tr [:blank:] "," | cut -d"," -f3)
  sed -i "s/$old_rate/$new_rate/" CastNet_parameters_${new_rate_name}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${new_rate_name}/" CastNet_${new_rate_name}.py
  date +%H,%M,%S > start_stop_${new_rate_name}-mutrate.csv
  python CastNet_${new_rate_name}.py ${new_rate_name}_indivsRun
  date +%H,%M,%S >> start_stop_${new_rate_name}-mutrate.csv
done
cd ~/Programs/CastNet
echo "Finished third bit successfully, DONE!"