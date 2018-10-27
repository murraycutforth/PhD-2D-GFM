
allres=( 10 20 40 80 160 320 640 )

for i in "${allres[@]}"
do
	cp ./../settingsfile.txt ./validation_settingsfiles/settingsfile_${i}.txt

	sed -i "s/Nx.*/Nx ${i}/g" ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i "s/Ny.*/Ny ${i}/g" ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i 's/FS.*/FS Godunov/g' ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i 's/test_case.*/test_case GDA/g' ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i "s/ITM.*/ITM LS5/g" ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/validation/#g' ./validation_settingsfiles/settingsfile_${i}.txt
	sed -i "s/GFM .*/GFM OGFM/g" ./validation_settingsfiles/settingsfile_${i}.txt
done

