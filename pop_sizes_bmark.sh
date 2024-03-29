#mkdir ~/Programs/CastNet/10-2K_PopSize_Bmarks
cd ~/Programs/CastNet/10-2K_PopSize_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
echo "pop_size,start,stop" > pop_sizes_times_parallel.csv
for new_pop_size in 10 50 100 200 300 400 500 1000 2000
do
  cp CastNet_parameters.py CastNet_parameters_${new_pop_size}.py
  cp CastNet.py CastNet_${new_pop_size}.py
  old_pop_size=$(egrep -o "pop_size = [0-9]+" CastNet_parameters_${new_pop_size}.py)
  sed -i "s/$old_pop_size/pop_size = $new_pop_size/" CastNet_parameters_${new_pop_size}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${new_pop_size}/" CastNet_${new_pop_size}.py
  start=$(date +%H:%M:%S)
  python CastNet_${new_pop_size}.py ${new_pop_size}_indivsRun
  stop=$(date +%H:%M:%S)
  echo $new_pop_size,$start,$stop >> pop_sizes_times_parallel.csv 
done
cd ~/Programs/CastNet
echo "Finished second bit successfully"
