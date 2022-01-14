#cp square_for_restart_od_17.8_ydn.1_ss90.nml input.nml
#./RUN1p
#cp ./RESTART/* ./S_17.8_ydn.1_ss90_r/

cp square_for_restart_od_17.8_ydn.1.nml input.nml
./RUN1p
cp ./RESTART/* ./S_17.8_ydn.1_r/
