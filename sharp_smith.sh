#!/bin/bash

conda activate smithHunter_env

#Default
organism="Unknown"
working_dir="$PWD"
mode="score"
path_bedfiles="NO"
path_smith="NO"
t_five_score=0.50
t_three_score=0.00
penalty=0.10
n_thre=0.50
sex="NA"



#options
while getopts "B:F:M:5:3:S:P:" opt; do
case "$opt" in
	B)path_bedfiles=$(echo $OPTARG | sed  "s/\/$//");;
	F)path_smith=$(echo $OPTARG | sed  "s/\/$//");;
	M)mode=$OPTARG;;
        5)t_five_score=$OPTARG;;
        3)t_three_score=$OPTARG;;
	S)sex=$OPTARG;;
	P)penalty=$OPTARG;;
        \?) echo -e "Argument Error in command line\n\nOptions:\n-M <mode: list-score-filter>\n-5 <five prime score>\n-3 <three prime score>\n-P <penalty>-B <path to bed files>\n-F <path to smith fasta file>\n-S <sex>"
       	exit 1;;
esac
done



if [[ $path_bedfiles == "NO" ]]; then
	echo "error in -B option. Folder of bed files must be specified";
	exit 1
fi

if [[ $mode == "filter" && $path_smith == "NO" ]]; then
        echo "error in -F option. In filter mode smith fasta file must be specified";
	exit 1 
fi


Rscript scripts/sharp_smith.R --mode=$mode --t_five_score=$t_five_score --t_three_score=$t_three_score --path_bedfiles=$path_bedfiles --path_smith=$path_smith --sex=$sex;

if [[ $mode == "filter" ]]; then
	mv $path_smith"/presumptive_smithRNAs.fa" $path_smith"/presumptive_smithRNAs_old.fa"
	cd $path_bedfiles;
	for i in $(cut -f 1 table.tpm); do 
		grep -A 1 "clusterid"$i"_" $path_smith"/presumptive_smithRNAs_old.fa" >> $path_smith"/presumptive_smithRNAs.fa"; 
	done

rm table.tpm;
cd $working_dir;

fi

conda deactivate 

exit

