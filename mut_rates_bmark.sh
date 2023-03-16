cd ~/Programs/CastNet/1e-8_1e-1_MutRate_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
echo "rate,start,stop" > mut_rates_times_parallel.csv
for line in {1..14}
do
  new_rate_name=$(sed -n ${line}p rates_table_full.list | cut -d"," -f1)
  new_rate=$(sed -n ${line}p rates_table_full.list | cut -d"," -f2)
  cp CastNet_parameters.py CastNet_parameters_${new_rate_name}.py
  cp CastNet.py CastNet_${new_rate_name}.py
  old_rate=$(grep seq_mutation_rate CastNet_parameters_${new_rate_name}.py | tr [:blank:] "," | cut -d"," -f3)
  sed -i "s/$old_rate/$new_rate/" CastNet_parameters_${new_rate_name}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${new_rate_name}/" CastNet_${new_rate_name}.py
  start=$(date +%H:%M:%S)
  python CastNet_${new_rate_name}.py ${new_rate_name}_indivsRun
  stop=$(date +%H:%M:%S)
  echo $new_rate,$start,$stop >> mut_rates_times_parallel.csv
done
cd ~/Programs/CastNet
echo "Finished mutation rate benchmarks successfully."
