
allGFM=( MGFM realGFM riemannGFM )

for i in "${allGFM[@]}"
do
	cp ./../settingsfile.txt ./usb_settingsfiles/settingsfile_200_${i}.txt

	sed -i 's/Nx.*/Nx 200/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/Ny.*/Ny 200/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/test_case.*/test_case underwater_shocked_bubble/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's/ITM.*/ITM CLSVOF/g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./usb_settingsfiles/settingsfile_200_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/underwater_shocked_bubble/#g' ./usb_settingsfiles/settingsfile_200_${i}.txt
	
	cp ./../settingsfile.txt ./usb_settingsfiles/settingsfile_400_${i}.txt

	sed -i 's/Nx.*/Nx 400/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/Ny.*/Ny 400/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/CFL.*/CFL 0.4/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/sim_type.*/sim_type twofluid/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/FS.*/FS MUSCL/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/test_case.*/test_case underwater_shocked_bubble/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's/ITM.*/ITM CLSVOF/g' ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i "s/GFM.*/GFM ${i}/g" ./usb_settingsfiles/settingsfile_400_${i}.txt
	sed -i 's#outputpath.*#outputpath /local/data2/public/mcc74/2D-GFM/underwater_shocked_bubble/#g' ./usb_settingsfiles/settingsfile_400_${i}.txt
done
