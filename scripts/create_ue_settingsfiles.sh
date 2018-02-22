
allGFM=( MGFM realGFM riemannGFM )

for i in "${allGFM[@]}"
do
	cp ./../settingsfile.txt ./ue_settingsfiles/settingsfile_400_${i}.txt

	sed -i 's/Nx.*/Nx 400/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/Ny.*/Ny 400/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/test_case.*/test_case underwater_explosion/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/ITM.*/ITM LS9/g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./ue_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/underwater_explosion/#g' ./ue_settingsfiles/settingsfile_400_${i}.txt
	
	cp ./../settingsfile.txt ./ue_settingsfiles/settingsfile_200_${i}.txt

	sed -i 's/Nx.*/Nx 200/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/Ny.*/Ny 200/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/test_case.*/test_case underwater_explosion/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/ITM.*/ITM LS9/g' ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./ue_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/underwater_explosion/#g' ./ue_settingsfiles/settingsfile_200_${i}.txt
done
