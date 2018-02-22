
allGFM=( OGFM MGFM realGFM riemannGFM )

for i in "${allGFM[@]}"
do
	cp ./../settingsfile.txt ./he_test_3_settingsfiles/settingsfile_325_${i}.txt

	sed -i 's/Nx.*/Nx 325/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/Ny.*/Ny 89/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/test_case.*/test_case shocked_helium_bubble/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's/ITM.*/ITM CLSVOF/g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/he_test_3/#g' ./he_test_3_settingsfiles/settingsfile_325_${i}.txt
	
	
	cp ./../settingsfile.txt ./he_test_3_settingsfiles/settingsfile_650_${i}.txt

	sed -i 's/Nx.*/Nx 650/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/Ny.*/Ny 178/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/test_case.*/test_case shocked_helium_bubble/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's/ITM.*/ITM CLSVOF/g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/he_test_3/#g' ./he_test_3_settingsfiles/settingsfile_650_${i}.txt
done
