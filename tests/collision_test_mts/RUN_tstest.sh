#!/bin/bash

rm iceberg_trajectories.nc

space="s_"
nc="c_.nc"
r="./results/iceberg_trajectories_"

s1='10' # 960' #'10 30 60 90'
c1='1.e-5' # 1.e-4'

for i in $s1
do
	for j in $c1
	do
	    	result="$r$i$space$j$nc"
		echo $i
		echo $j
		echo $result

		sed "s/<sec>/$i/; s/<coef>/$j/g" \
		    input_tstest.nml > input.nml

		srun -n4 ./../../build/bergs.x

		mv iceberg_trajectories.nc $result
	done
done





