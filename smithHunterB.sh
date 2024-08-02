#######################################################################################################################
#                                                                                                                     #
# SMITHHUNTER is designed to identify putative smithRNAs in a species of interest starting from small RNA libraries,  #
# the sequence of the mitochondrial genome, the sequence of the transcriptome inclusive of UTR annotations and,       #
# optionally, the sequence of the nuclear genome. The first module focuses on the identification and filtering of     #
# presumptive smithRNA sequences, defined as centroids of clusters with significant transcription levels and a narrow #
# 5’ transcription boundary. The second module deals with the identification of possible nuclear targets and          #
# pre-miRNA-like precursor structures for presumptive smithRNAs. A third script is provided to help identify          #
# smithRNAs with narrow start/endpoints.                                                                              #
#                                                                                                                     #
# Copyright (C) 2024 Giovanni Marturano, Diego Carli.                                                                 #
#                                                                                                                     #
# This program is free software: you can redistribute it and/or modify                                                #
# it under the terms of the GNU General Public License as published by                                                #
# the Free Software Foundation, either version 3 of the License, or	                                              #
# (at your option) any later version.                                                                                 #
#                                                                                                                     #
# This program is distributed in the hope that it will be useful,                                                     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                       #
# GNU General Public License for more details.                                                                        #
#                                                                                                                     #
# You should have received a copy of the GNU General Public License                                                   #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                               #
#                                                                                                                     #
#######################################################################################################################

#!/bin/bash

conda activate smithHunter_env

# Calculate terminal width
term_width=$(tput cols)
text_width=68

# Calculate left padding for centering
left_padding=$(( (term_width - text_width) / 2 ))


# Print centered text
printf "%*s${GREEN}\033[1m                    _ __  __    __  ____  ___   __________________ \n" $left_padding
printf "%*s${GREEN}\033[1m    _________ ___  (_) /_/ /_  / / / / / / / | / /_  __/ ____/ __ \ \n" $left_padding
printf "%*s${GREEN}\033[1m   / ___/ __ '__ \/ / __/ __ \/ /_/ / / / /  |/ / / / / __/ / /_/ / \n" $left_padding
printf "%*s${GREEN}\033[1m  (__  ) / / / / / / /_/ / / / __  / /_/ / /|  / / / / /___/ _, _/  \n" $left_padding
printf "%*s${GREEN}\033[1m /____/_/ /_/ /_/_/\__/_/ /_/_/ /_/\____/_/ |_/ /_/ /_____/_/ |_|  \n${NC}" $left_padding

########### OPTIONS #########################################
#Default
home=$PWD
organism="Unknown"
seed_start="4"
seed_end="10"
seed_mismatch="0"
RNAfoldT=25
dist1=15
dist2=50

while getopts "W:O:P:X:Y:m:R:1:2:" opt; do
case "$opt" in
	W)home=$(echo $OPTARG | sed  "s/\/$//");;
        O)organism=$OPTARG;;
        P)PITA=$(echo $OPTARG | sed  "s/\/$//");;
        X)seed_start=$OPTARG;;
        Y)seed_end=$OPTARG;;
	m)seed_mismatch=$OPTARG;;
	R)RNAfoldT=$OPTARG;;
	1)dist1=$OPTARG;;
        2)dist2=$OPTARG;;

        \?) echo -e "Argument Error in command line\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1;;
esac
done


im_here=$(pwd)

#change relative in absolute paths
if [[ -d "$home" ]]; then
        if [[ "$home" != /* ]]; then
                cd $home;
                home=$(pwd)
                cd $im_here
        fi
else
        echo ""
        echo "Working directory not found"
        exit 1

fi

if [[ -d "$PITA" ]]; then
        if [[ "$PITA" != /* ]]; then
                cd $PITA;
                PITA=$(pwd)
                cd $im_here
        fi

else
        echo ""
        echo "PITA folder not found"
        exit 1

fi





smithRNA_fasta=$home/"7_"$organism"_smithRNAs"
BLAST=$home/"8_"$organism"_BLAST"
thermodynamics=$home/"9_"$organism"_thermodynamics"
targets=$home/"10_"$organism"_targets"
folding=$home/"11_"$organism"_folding"
out_mit=$home/"0_"$organism"_mit"/$organism"_mit.fasta"
main_outputs=$home/$organism"_main_outputs"


Transcript=$home/$organism"_Transcripts.fasta"
UTR_bed=$home/$organism"_UTR.bed"

#controllo gli argomenti
echo ""
echo "Working directory=$home"
echo ""
echo "Organims=$organism"
echo ""

if [ -f $smithRNA_fasta/presumptive_smithRNAs.fa ] ; then
        echo "presumptive smithRNAs fasta file OK"
        echo ""
else
        echo -e "presumptive smithRNAs fasta file missing. It should be located in the $smithRNA_fasta directory and named presumptive_smithRNAs.fa\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi

if [ -f $out_mit ] ; then
        echo "MITO genome OK"
        echo ""
else
        echo -e "MITO-genome fasta file missing. It should be located in the $home"/0_"$organism"_mit" directory and named $organism"_mit.fasta"\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi


if [ -f $Transcript ] ; then
        echo "Transcripts fasta file OK"
        echo ""
else
        echo -e "Transcripts fasta file missing. It should be located in the working directory and named $organism"_Transcripts.fasta"\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi


if [ -f $UTR_bed ] ; then
        echo "UTR bed file OK"
        echo ""
else
        echo -e "UTR bed file missing. It should be located in the working directory and named $organism"_UTR.bed" \n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi

if [ -f $PITA/pita_prediction.pl ] ; then
        echo "PITA OK"
        echo ""
else
        echo -e "PITA software missing\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi

#List of Seq IDs in the fasta transcripts file
grep ">" $Transcript | sed 's/>//g' > FAseqIDs

#List of Seq IDs in the UTR bed file
cut -f 1 $UTR_bed > BEDseqIDs

#Number of Seq IDs in UTR file
N=$(cat BEDseqIDs | sort | uniq | wc -l)

#Number of identical Seq IDs in the bed and fasta file..must be equal to N
N2=$(grep -w -f BEDseqIDs FAseqIDs | wc -l)

if [ "$N" == "$N2" ]; then
        echo "sequences IDs in $Transcript  and $UTR_bed files OK"
        echo ""
else
        echo "sequences IDs in $Transcript and $UTR_bed files do not correspond"
        echo ""
        exit 1
fi

rm *seqIDs

if [[ $seed_start -gt 0 && $seed_start -le 20 ]] ; then
        echo "seed start=$seed_start"
        echo ""
else
        echo -e "Error in -X opion. Values between 1 and 20 are accepted\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi


if [[ $seed_end -gt 0 && $seed_end -le 20 && $seed_end -gt $seed_start ]] ; then
        echo "seed end=$seed_end"
        echo ""
else
        echo -e "Error in -Y opion. Values between 1 and 20 are accepted. -Y must be higher than -X\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi

if (( $(echo "$seed_mismatch >= 0 && $seed_mismatch <= 2" | bc -l) )); then
        echo "seed mismatches allowed=$seed_mismatch"
        echo ""
else
        echo -e "Error in -m opion. From 0 to 2 miss-matches allowed in the seed sequence\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
        exit 1
fi

if [[ $dist1 -gt 0 && $dist1 -le 100 ]]; then
	echo "Pre-smithRNA distance 1=$dist1"
	echo ""
else
	echo -e "Error in -1 opion. Values between 1 and 100 are accepted\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
	exit 1
fi

if [[ $dist2 -gt 0 && $dist2 -le 100 ]]; then
        echo "Pre-smithRNA distance 2=$dist2"
	echo ""
else
        echo -e "Error in -2 opion. Values between 1 and 100 are accepted\n\nOptions:\n-W <working dir path>\n-P <pita software path>\n-O <organism ID>\n-X <start seed position>\n-Y <end seed positon>\n-m <mismatches allowed into seed region>\n-R <RNA folding temperature>\n-1 <pre-smithRNA distance1>\n-2 <pre-smithRNA distance2>"
	exit 1
fi

echo "RNA folding temperature=$RNAfoldT"
echo ""
sleep 5

##########FOLDERS CREATION #############################

if [ -f $BLAST ]
        then rm -rf $BLAST
        else mkdir $BLAST
fi

if [ -f $thermodynamics ]
        then rm -rf $thermodynamics
        else mkdir $thermodynamics
fi

if [ -f $targets ]
        then rm -rf $targets
        else mkdir $targets
fi

if [ -f $folding ]
        then rm -rf $folding
	else mkdir $folding
fi

#UTR extraction
bedtools getfasta -fi $Transcript -bed $UTR_bed  > $home/$organism"_UTR.tmp.fasta"


#formatting UTRs file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $home/$organism"_UTR.tmp.fasta" > $home/$organism"_UTR.fasta"

rm $home/$organism"_UTR.tmp.fasta"

#reversing seed and whole smithRNAs
revseq $smithRNA_fasta/presumptive_smithRNAs.fa -tag FALSE -sbegin1 $seed_start -send1 $seed_end -outseq $smithRNA_fasta/RC_seed_smithRNAs.fa
revseq $smithRNA_fasta/presumptive_smithRNAs.fa -tag FALSE -outseq $smithRNA_fasta/RC_smithRNAs.fa

cd $BLAST

if [ -f BLASTDBs ]
        then rm BLASTDBs
        else mkdir BLASTDBs
fi

if [ -f Whole_UTRs ]
        then rm Whole_UTRs
        else mkdir Whole_UTRs
fi

if [ -f Selected_UTRs ]
        then rm Selected_UTRs
        else mkdir Selected_UTRs
fi

if [ -f Queries ]
        then rm Queries
        else mkdir Queries
fi

#work indipendently on each smithRNA
num_smith=$(grep -c ">" $smithRNA_fasta/presumptive_smithRNAs.fa)

echo -e "SMITH\tTARGET\tPITA_dG\tRNAhybrid_dG" > $thermodynamics/dG_INFO.tmp.txt
echo -e "SMITH\tTARGET\tPITA_dG\tRNAhybrid_dG" > $thermodynamics/dG_INFO.PASS.txt

for s in $(seq 1 $num_smith)
        do curr_smith=$(head -n $(( s*2-1 )) $smithRNA_fasta/presumptive_smithRNAs.fa | tail -n 1 | sed 's/>//g')
        if [ -f $curr_smith.BLAST.target ]
                then rm $curr_smith.BLAST.target
                fi
        if [ -f $curr_smith.selected.UTR.fa ]
                then rm $curr_smith.selected.UTR.fa
                fi
        if [ -f $curr_smith.targets ]
                then rm $curr_smith.targets
                fi
	
	curr_seed=$(head -n $(( s*2 )) $smithRNA_fasta/RC_seed_smithRNAs.fa | tail -n 1)
	curr_RNA=$(head -n $(( s*2 )) $smithRNA_fasta/RC_smithRNAs.fa | tail -n 1)

	#grepping seed in the UTR file, and edit it till a regular fasta file
        grep_sequences=$(agrep -$seed_mismatch "$curr_seed" $home/$organism"_UTR.fasta");
        for i in $grep_sequences; do
		grep  -w -B 1 "$i" $home/$organism"_UTR.fasta" >> Whole_UTRs/$curr_smith".seed.tmp";
        done

	#removing eventual duplicates from target file
	awk '/^>/{f=!d[$1];d[$1]=1}f' Whole_UTRs/$curr_smith".seed.tmp" > Whole_UTRs/$curr_smith".seed";

	rm  Whole_UTRs/$curr_smith".seed.tmp";


        all_targets=$(grep ">" Whole_UTRs/$curr_smith".seed" | sed 's/>//g')


	#build a blastn database on the seed file (with all the UR with a match with the curr_seed, and keep only if there is a match of at least 11 nucleotides.
        makeblastdb -in Whole_UTRs/$curr_smith".seed" -dbtype nucl -title $curr_smith"_UTR" -out BLASTDBs/$curr_smith"_UTR"
        head -n $(( s*2 )) $smithRNA_fasta/presumptive_smithRNAs.fa | tail -n 2 > Queries/$curr_smith".query.fa"
        blastn -task blastn-short -db BLASTDBs/$curr_smith"_UTR" -query Queries/$curr_smith".query.fa" -strand minus -evalue 10000  -max_target_seqs 100000 -out $curr_smith".temp.BLAST.target" -outfmt '6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send sstrand evalue'

        num_lines=$(wc -l $curr_smith".temp.BLAST.target" | awk '{print $1}')
        for h in `seq 1 $num_lines`
                do if [ $(head -n $h $curr_smith".temp.BLAST.target" | tail -n 1 | awk '{print $4}') -ge 11 ]
                        then head -n $h $curr_smith".temp.BLAST.target" | tail -n 1 >> $curr_smith".BLAST.target"
                        fi
                done

        #retrieve every UTRs that has passed the prevoius two tests (complete alignment of seed and 11 match on the whole smithRNA)
	if [ -f $curr_smith'.BLAST.target' ] ; then
        for U in `awk '{print $2}' $curr_smith'.BLAST.target'`
                do grep -A 1 -w $U $home/$organism"_UTR.fasta" >> Selected_UTRs/$curr_smith'.selected.UTR.tmp.fa'
                sed -i 's/\([^\t]*\).*/\1/g' Selected_UTRs/$curr_smith'.selected.UTR.tmp.fa'
                cat Selected_UTRs/$curr_smith'.selected.UTR.tmp.fa' | paste - - | sort -u | awk '{print $1"\n"$2}' > Selected_UTRs/$curr_smith'.selected.UTR.fa'
                selected_targets=$(grep ">" Selected_UTRs/$curr_smith'.selected.UTR.fa' | sed 's/>//g')
        done

	rm Selected_UTRs/$curr_smith'.selected.UTR.tmp.fa'

	fi

        #running PITA and RNAhybrid
	cd $PITA
        perl pita_prediction.pl -utr $BLAST/Selected_UTRs/$curr_smith'.selected.UTR.fa' -mir $BLAST/Queries/$curr_smith'.query.fa' -flank_up 3 -flank_down 15 -prefix $thermodynamics/$curr_smith
        if [ -s $BLAST/Selected_UTRs/$curr_smith'.selected.UTR.fa' ]
                then
#               #n option is the maximum allowed length of the query file (default 30), m option the same, but for the target (default 2000)
#               #s flag specify the target dataset to assume as for the quick estimate of extreme value distribution parameters. Available option are 3utr_fly, 3utr_worm and 3utr_human
                RNAhybrid -f 3,10 -e -15 -n 100 -m 1000000 -s 3utr_fly -t $BLAST/Selected_UTRs/$curr_smith'.selected.UTR.fa' -q $BLAST/Queries/$curr_smith'.query.fa' > $thermodynamics/$curr_smith'.RNAhybrid_results.out'
        fi

        for target in $selected_targets; do
                DG_Pita=$(grep --no-messages -w "$target" $thermodynamics/$curr_smith"_pita_results.tab" | cut -f 7 | head -n 1)
                DG_RNAhybrid=$(grep --no-messages -A 6 -w "$target" $thermodynamics/$curr_smith".RNAhybrid_results.out" | grep "mfe" | sed 's/mfe: //g' | sed 's/ /\t/g' | cut -f 1 | head -n 1)
                echo -e "$curr_smith\t$target\t$DG_Pita\t$DG_RNAhybrid" >> $thermodynamics/dG_INFO.tmp.txt
        done

	cd $BLAST

done 


cd $thermodynamics

awk '$3 != "" || $4 != "" {print}' dG_INFO.tmp.txt > dG_INFO.txt
rm dG_INFO.tmp.txt

awk '{OFS="\t"; $3=int($3); $4=int($4); print}' dG_INFO.txt > dG_INFO.tmp2.txt

awk '$3 <= -9 && $4 <= -15 {print}' dG_INFO.tmp2.txt >> dG_INFO.PASS.txt

rm dG_INFO.tmp2.txt


selected_smiths=$(cut -f 1 dG_INFO.PASS.txt | sort | uniq | grep -v "SMITH")

for smith in $selected_smiths; do
       selected_targets=$(grep -w "$smith" dG_INFO.PASS.txt | cut -f 2 | sed s'/:/\t/g' | cut -f 1)
        for target in $selected_targets; do
                grep -w -A 1 "$target" $Transcript >> $targets/$smith".targets"
		#removing duplicates from target file
                cat $targets/$smith'.targets' | paste - - > $targets/tab_$smith'.targets'
                sort -u $targets/tab_$smith'.targets' > $targets/'sorted_tab_'$smith'.targets'
                awk '{print $1"\n"$2}' $targets/'sorted_tab_'$smith'.targets' > $targets/$smith'.target.fa'

        done
done


#adding headers to BLAST.target files
sed -i '1 i qseqid\tsseqid\tpident\tnident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue' $BLAST/*.BLAST.target

cd $targets

#final renaming and stuff
if [ -f validated.smithRNAs.fa ]
        then rm validated.smithRNAs.fa
        fi


for validated in $(ls  *.target.fa | sed 's/\.target\.fa//g'); do
	grep -A 1 -w "$validated" $smithRNA_fasta/presumptive_smithRNAs.fa >> $smithRNA_fasta/candidate_smithRNAs.fa;
done


#removing temporary files
rm -r $BLAST/*.temp.BLAST.target $BLAST/BLASTDBs $BLAST/Queries $BLAST/Selected_UTRs  $BLAST/Whole_UTRs
rm $targets/*targets

##########FOLDING#################################################################

grep ">" $smithRNA_fasta/candidate_smithRNAs.fa | sed 's/_/\t/g' | sed 's/size//g' | sed 's/strand//g' | sed 's/pos//g' > $folding/list.txt

cd $folding

mkdir plots_folding

echo -e "Pre_smith_ID\tdG1\tdG2" > dG.txt


num_pre_smith=$(cat list.txt | wc -l)
pre_smith=1

while  [ $pre_smith -le $num_pre_smith ]; do
        cluster=$(head -n $pre_smith list.txt | cut -f 1 | tail -n 1 | sed 's/>//g')
        size=$(head -n $pre_smith list.txt | cut -f 2 | tail -n 1)
        start=$(head -n $pre_smith list.txt | cut -f 3 | tail -n 1)
        end=$(head -n $pre_smith list.txt | cut -f 4 | tail -n 1)
        strand=$(head -n $pre_smith list.txt | cut -f 5 | tail -n 1)
        ID=$(echo $cluster"_"$size"_"$start"_"$end"_"$strand)

        if [ $strand == "-" ];
                then
                        start1=$(($start - $dist1))
                        end1=$(($end + $dist2))
                        start2=$(($start - $dist2))
                        end2=$(($end + $dist1))
                        seqret -sequence $out_mit -sbegin $start1 -send $end1 -outseq $ID".pre1.fa"
                        seqret -sequence $out_mit -sbegin $start2 -send $end2 -outseq $ID".pre2.fa"
                        revseq $ID".pre1.fa" -tag FALSE -outseq "RC_"$ID".pre1.fa"
                        revseq $ID".pre2.fa" -tag FALSE -outseq "RC_"$ID".pre2.fa"
                        RNAfold -T $RNAfoldT -i "RC_"$ID".pre1.fa" > "RC_"$ID".pre1.dG"
                        RNAfold -T $RNAfoldT -i "RC_"$ID".pre2.fa" > "RC_"$ID".pre2.dG"

                        sed -i "s/>.*/>1$ID/g" $ID".pre1.fa"
                        sed -i "s/>.*/>2$ID/g" $ID".pre2.fa"
                        sed -i "s/>.*/>1$ID/g" "RC_"$ID".pre1.fa"
                        sed -i "s/>.*/>2$ID/g" "RC_"$ID".pre2.fa"
                        sed -i "s/>.*/>1$ID/g" "RC_"$ID".pre1.dG"
                        sed -i "s/>.*/>2$ID/g" "RC_"$ID".pre2.dG"

                        cat "RC_"$ID".pre1.dG" | RNAplot -o svg
                        cat "RC_"$ID".pre2.dG" | RNAplot -o svg


                        DG1=$(sed -e 's/  */\n/g' "RC_"$ID".pre1.dG" | tail -n 1 | sed s'/(//g' | sed s'/)//g' | sed 's/\..*//g')
                        DG2=$(sed -e 's/  */\n/g' "RC_"$ID".pre2.dG" | tail -n 1 | sed s'/(//g' | sed s'/)//g' | sed 's/\..*//g')

                        echo -e "$ID\t$DG1\t$DG2" >> dG.txt


                else
                        start1=$(($start - $dist1))
                        end1=$(($end + $dist2))
                        start2=$(($start - $dist2))
                        end2=$(($end + $dist1))
                        seqret -sequence $out_mit -sbegin $start1 -send $end1 -outseq $ID".pre1.fa"
                        seqret -sequence $out_mit -sbegin $start2 -send $end2 -outseq $ID".pre2.fa"
                        RNAfold -T $RNAfoldT -i $ID".pre1.fa" > $ID".pre1.dG"
                        RNAfold -T $RNAfoldT -i $ID".pre2.fa" > $ID".pre2.dG"

                        sed -i "s/>.*/>1$ID/g" $ID".pre1.fa"
                        sed -i "s/>.*/>2$ID/g" $ID".pre2.fa"
                        sed -i "s/>.*/>1$ID/g" $ID".pre1.dG"
                        sed -i "s/>.*/>2$ID/g" $ID".pre2.dG"

                        cat $ID".pre1.dG" | RNAplot -o svg
                        cat $ID".pre2.dG" | RNAplot -o svg

                        DG1=$(sed -e 's/  */\n/g' $ID".pre1.dG" | tail -n 1 | sed s'/(//g' | sed s'/)//g' | sed 's/\..*//g')
                        DG2=$(sed -e 's/  */\n/g' $ID".pre2.dG" | tail -n 1 | sed s'/(//g' | sed s'/)//g' | sed 's/\..*//g')


                        echo -e "$ID\t$DG1\t$DG2" >> dG.txt


        fi

        pre_smith=$(( pre_smith + 1 ))
done

mv *svg plots_folding

rm $folding/list.txt


#copy the main outputs in the fìdirectory
mkdir $main_outputs/targets
cp $smithRNA_fasta/candidate_smithRNAs.fa $main_outputs
cp $thermodynamics/dG_INFO.PASS.txt $main_outputs
cp -r $folding/plots_folding $main_outputs/Plots 
cp $targets/*fa $main_outputs/targets


conda deactivate 

exit 0

