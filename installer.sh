#!/bin/bash

echo -e "SmithHunter installer.This script allows you to create smithHunter conda enviroment and install PITA software ...\n"
sleep 3

#find conda
conda=$(which conda | wc -l)
pita=$(find /home/$(whoami) -name pita_prediction.pl | wc -l)

#make SmithHunter executable:
chmod 755 smithHunter.yml
chmod 755 smithHunterA.sh
chmod 755 smithHunterB.sh
chmod 755 scripts/esplora2.R
chmod 755 scripts/Makeplots.R

#check conda 
if [[ $conda == 0 ]]; then
	echo -e "Conda program not found. Please install conda to proceed\n"
	echo -e "Do you want us to install Conda for you? [y/n]\n"
	read yesorno
	if [[ $yesorno == "n" ]]; then
		echo -e "\nPlease install conda to proceed"
		exit 1
	elif [[ $yesorno == "y" ]]; then
		echo -e "\nConda installation...\n"
		# download the installer
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P /home/$(whoami)
		# make it executable
		chmod 755 /home/$(whoami)/Miniconda3-latest-Linux-x86_64.sh
		# run the installer
		/home/$(whoami)/Miniconda3-latest-Linux-x86_64.sh
		echo -e "Creation of smithHunter conda enviroment\n"
		sleep 3
		# create conda environment
		conda env create -f smithHunter.yml
		condaenv=$(find /home/$(whoami) -name conda | grep "smithHunter_env" | wc -l)
		if [[ $condaenv -gt 0 ]]; then
			echo -e "\ndone\n"
		else 
			echo "error in env creation"
			exit 1
		fi
	else
		echo -e "\nonly "y" or "n" are accepted as answer";
		exit 1
	fi
elif [[ $conda -gt 0 ]]; then
        echo -e "Conda program found\n\nCreation of smithHunter conda enviroment\n"
        sleep 3
	# create conda environment
	conda env create -f smithHunter.yml
	condaenv=$(find /home/$(whoami) -name conda | grep "smithHunter_env" | wc -l)
	if [[ $condaenv -gt 0 ]]; then
		echo -e "\ndone\n"
	else
		echo "error in env creation"
		exit 1
        fi
fi

#check pita
if [[ $pita == 0 ]]; then
        echo -e "\nPITA program not found. Please install PITA to proceed\n"
        echo "Do you want us to install PITA for you? [y/n]"
        read yesorno2
        if [[ $yesorno2 == "n" ]]; then
                echo -e "\nPlease install PITA to proceed"
                exit 1
        elif [[ $yesorno2 == "y" ]]; then
                echo -e "\nPITA installation...\n"
                # download and install PITA
		git clone https://github.com/elenacorni/pita
		cd pita/pita_cpp/
		make install
		if [[ -f pita_prediction.pl ]]; then
			echo -e "\nPITA software installed in $PWD"
			sleep 2
			echo -e "\ninstallation complete"
		fi
        else
                echo -e "\nonly "y" or "n" are accepted as answer";
        fi
elif [[ $pita -gt 1 ]]; then
        echo -e "PITA program found\n"
	sleep 2
	echo -e "\ninstallation complete"
fi;


