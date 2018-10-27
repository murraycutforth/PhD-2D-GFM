
allLS=( LS1 LS5 CLSVOF EMOF2 )

for i in "${allLS[@]}"
do
	cp ./../settingsfile.txt ./tstm_settingsfiles/settingsfile_${i}.txt

	sed -i 's/Nx.*/Nx 210/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's/Ny.*/Ny 90/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's/test_case.*/test_case TSTM/g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i "s/ITM.*/ITM ${i}/g" ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/tstm/#g' ./tstm_settingsfiles/settingsfile_${i}.txt
	sed -i "s/GFM .*/GFM OGFM/g" ./tstm_settingsfiles/settingsfile_${i}.txt
done

sed -i "s/GFM OGFM/GFM OGFMVOF/g" ./tstm_settingsfiles/settingsfile_EMOF2.txt
