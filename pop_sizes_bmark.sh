cd ~/Programs/CastNet/10-1K_PopSize_Bmarks
cp ../CastNet.py ../CastNet_parameters.py ../CastNet_out_funcs.py .
for new_pop_size in 10 50 100 200 300 400 500 1000 2000
do
  cp CastNet_parameters.py CastNet_parameters_${new_pop_size}.py
  cp CastNet.py CastNet_${new_pop_size}.py
  old_pop_size=$(egrep -o "pop_size = [0-9]+" CastNet_parameters_${new_pop_size}.py)
  sed -i "s/$old_pop_size/pop_size = $new_pop_size/" CastNet_parameters_${new_pop_size}.py
  sed -i "s/CastNet_parameters/CastNet_parameters_${new_pop_size}/" CastNet_${new_pop_size}.py
  date +%H,%M,%S > start_stop_${new_pop_size}-indivs.csv
  python CastNet_${new_pop_size}.py ${new_pop_size}_indivsRun
  date +%H,%M,%S >> start_stop_${new_pop_size}-indivs.csv
done
cd ~/Programs/CastNet
echo "Finished second bit successfully"
