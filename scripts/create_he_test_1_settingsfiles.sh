
allLS=( LS1 LS5 LS9 CLSVOF )

for i in "${allLS[@]}"
do
	cp ./../settingsfile.txt ./he_test_1_settingsfiles/settingsfile_${i}.txt

	sed -i 's/Nx.*/Nx 325/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's/Ny.*/Ny 89/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's/test_case.*/test_case shocked_helium_bubble/g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i "s/ITM.*/ITM ${i}/g" ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i 's#outputpath.*#outputpath ./output/he_test_1/noreinit/#g' ./he_test_1_settingsfiles/settingsfile_${i}.txt
	sed -i "s/GFM.*/GFM MGFM/g" ./he_test_1_settingsfiles/settingsfile_${i}.txt
done
