
allLS=( LS1 LS5 CLSVOF EMOF2 )

for i in "${allLS[@]}"
do
	cp ./../settingsfile.txt ./shocked_r22_settingsfiles/settingsfile_${i}.txt

	sed -i 's/Nx.*/Nx 3560/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's/Ny.*/Ny 712/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's/test_case.*/test_case shocked_R22_bubble/g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i "s/ITM.*/ITM ${i}/g" ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/shocked_R22/#g' ./shocked_r22_settingsfiles/settingsfile_${i}.txt
	sed -i "s/GFM .*/GFM OGFM/g" ./shocked_r22_settingsfiles/settingsfile_${i}.txt
done

sed -i "s/GFM OGFM/GFM OGFMVOF/g" ./shocked_r22_settingsfiles/settingsfile_EMOF2.txt
