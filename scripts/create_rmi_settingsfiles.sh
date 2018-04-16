
allLS=( LS1 LS5 CLSVOF EMOF2 )

for i in "${allLS[@]}"
do
	cp ./../settingsfile.txt ./rmi_settingsfiles/settingsfile_${i}.txt

	sed -i 's/Nx.*/Nx 2400/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's/Ny.*/Ny 600/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's/test_case.*/test_case RMI/g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i "s/ITM.*/ITM ${i}/g" ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/rmi/#g' ./rmi_settingsfiles/settingsfile_${i}.txt
	sed -i "s/GFM .*/GFM OGFM/g" ./rmi_settingsfiles/settingsfile_${i}.txt
done

sed -i "s/GFM OGFM/GFM OGFMVOF/g" ./rmi_settingsfiles/settingsfile_EMOF2.txt
