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