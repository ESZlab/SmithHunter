#!/bin/bash

#need to copy manually in the working directory the nuclear and mitochondrial genomes of interest named <organism_nuc.fa> and <organism_mit.fasta> respectively
#need to run with bash -i to fix conda activate bias


########### OPTIONS #########################################
#Default
home=$PWD
fastq_folder=$home/fastq
scripts=scripts
organism="Unknown"
Trimming=SE
Shift=1
clusteringID=0.95
stringency=0.50
threads=5
min_rep=1
adapter_sequence=NA;
adapter_sequenceR2=NA;
script1=$(find /home/$(whoami) -name esplora2.R | grep "SmithHunter" | head -n 1)
Nscript1=$(find /home/$(whoami) -name esplora2.R | grep "SmithHunter" | wc -l)
script2=$(find /home/$(whoami) -name Makeplots.py | grep "SmithHunter" | head -n 1)
Nscript2=$(find /home/$(whoami) -name Makeplots.py | grep "SmithHunter" | wc -l )

#check R and Python scripts
if [[ $Nscript1 == 0 ]]; then
	echo ""
        echo "esplora2.R script missing, it must be present in the SmithHunter_v*/scripts folder"
        exit 1
elif [[ $Nscript1 -gt 1 ]]; then
	echo ""
        echo "esplora2.R present in more folders. This script will be used in the current analyses: $script1"
	sleep 4
fi

if [[ $Nscript2 == 0 ]]; then
	echo ""
        echo "Makeplots.py script missing, it must be present in the SmithHunter_v*/scripts folder"
        exit 1
elif [[ $Nscript2 -gt 1 ]]; then
	echo ""
        echo "Makeplots.py present in more folders. This script will be used in the current analyses: $script2"
	sleep 4
fi


while getopts "W:F:O:T:I:S:t:M:a:A:s:" opt; do
case "$opt" in
        W)home=$(echo $OPTARG | sed  "s/\/$//");;
        F)fastq_folder=$(echo $OPTARG | sed  "s/\/$//");;
        O)organism=$OPTARG;;
        T)Trimming=$OPTARG;;
        I)clusteringID=$OPTARG;;
        S)stringency=$OPTARG;;
        t)threads=$OPTARG;;
	M)min_rep=$OPTARG;;
	a)adapter_sequence=$OPTARG;;
	A)adapter_sequenceR2=$OPTARG;;
	s)Shift=$OPTARG;;
        \?) echo -e "Argument Error in command line\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1;;
esac
done

#Options Check
echo ""
echo "Working directory=$home"
echo ""


if ls "$fastq_folder"/*_[1-2].fastq.gz 1> /dev/null 2>&1 ; then
        echo "Raw fastq files directory=$fastq_folder"
        echo ""
else
        echo -e "fastq files missing. SmallRNA sequence file names should end in ‘_1.fastq.gz’ and ‘_2.fastq.gz’ for PE data and in ‘_1.fastq.gz’ for SE data\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T <trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi

if [ -f $home/$organism"_mit.fasta" ]; then
        echo ""Organism=$organism""
        echo ""
else
        echo -e "MITO fasta genome missing. The sequence of the mito-genome should be located in the working directory and named $organism"_mit.fasta"\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi


if [[ "$Trimming" == "PE" ]] ; then
	if ls "$fastq_folder"/*2.fastq.gz 1> /dev/null 2>&1 ; then
	        echo "Trimming in PE mode"
		echo ""
	else 
		echo "R2 fastq files missing"
		exit 1
	fi
elif  
	[[ "$Trimming" == "PE" && "$adapter_sequence" == "NA" && "$adapter_sequenceR2" != "NA" ||  "$Trimming" == "PE" && "$adapter_sequence" != "NA" && "$adapter_sequenceR2" == "NA" ]] ; then
	echo -e "In PE mode both or neither adapters sequences have to be specified\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1

elif 
	[[ "$Trimming" == "SE" && "$adapter_sequence" != "NA" ]] ; then
        echo "Trimming in SE mode"
        echo ""
elif 
	[[ "$Trimming" == "SE" && "$adapter_sequence" == "NA" ]] ; then
	echo -e "In SE mode adapters sequences have to be specified\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
	exit 1

elif 
        [[ "$Trimming" == "NO" ]] ; then
        echo "Trimming deactivated"
	echo ""

else
        echo -e "Error in -T arguments\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi

if (( $(echo "$clusteringID >= 0 && $clusteringID <= 1" | bc -l) )); then
        echo "Clustering ID=$clusteringID"
        echo ""
else
        echo -e "Error in -I opion. Accepted values betweem 0 and 1\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi

if (( $(echo "$stringency >= 0 && $stringency <= 1" | bc -l) )); then
        echo "Stringency=$stringency"
        echo ""
else
        echo -e "Error in -S opion. Accepted values between 0 and 1\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi


if [[ $threads -gt 0 && $threads -le 40 ]]; then
        echo "threads=$threads"
        echo ""
else
        echo -e "Error in -t opion. Accepted values between 1 and 40\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
        exit 1
fi

num_samples=$(ls $fastq_folder/*fastq.gz | sed s"/\_[1|2].fastq.gz//g" | sort | uniq | wc -l);

if [[ $min_rep -le $num_samples && $min_rep -ge 1 ]]; then
	echo "Min rep=$min_rep"
	echo ""
else
	echo -e "Error in -M opion. Accepted values between 1 and $num_samples\n\nOptions:\n-W <working dir path>\n-F <fastq dir path>\n-O <organism ID>\n-T<trimming:PE-SE-NO>\n-a <adapter sequence R1>\n-A <adapter sequence R2>\n-I <clustering identity:0-1>\n-S <stringency:0-1>\n-M <min replicates>\n-s <shift:mito-genome origin>\n-t <threads:1-40>"
	exit 1
fi

if [ -f $home/$organism"_nuc.fasta" ]; then
        echo "NUCLEAR fasta genome found"
        echo ""
else
        echo "NUCLEAR fasta genome NOT found, nuclear-mapping step will be skipped"
        echo ""
fi

#Folder and files definition
genome_mit=$home/"0_"$organism"_mit"
genome_nuc=$home/"0_"$organism"_nuc"
trimming=$home/"1_"$organism"_trimmed"
fastqc=$home/"2_"$organism"_fastqc"
alignments=$home/"3_"$organism"_alignments"
BedFiles=$home/"4_"$organism"_COV_and_BedFiles"
Clustering=$home/"5_"$organism"_clustering"
Plots=$home/"6_"$organism"_plots"
smithRNA_fasta=$home/"7_"$organism"_smithRNAs"
main_outputs=$home/$organism"_main_outputs"

#conda enviroment activation
conda activate smithHunterA_env

#folder creation
if [ -d $genome_mit ]
	then rm -rf $genome_mit
fi
mkdir $genome_mit

if [ -f $home/$organism"_nuc.fasta" ]; then
	if [ -d $genome_nuc ]
		then rm -rf $genome_nuc
	fi
	mkdir $genome_nuc
fi

if [ -d $trimming ]
	then rm -rf $trimming
fi

mkdir $trimming

if [ -d $fastqc ]
	then rm -rf $fastqc
fi

mkdir $fastqc

if [ -d $alignments ]
	then rm -rf $alignments
fi

mkdir $alignments

if [ -d $BedFiles ]
	then rm -rf $BedFiles
fi

mkdir $BedFiles

if [ -d $Clustering ]
	then rm -rf $Clustering
fi
mkdir $Clustering

if [ -d $Plots ]
	then rm -rf $Plots
fi

if [ -f $organism"_samples.txt" ]
        then rm $organism"_samples.txt"
fi

mkdir $Plots

if [ -f $smithRNA_fasta ]
        then rm $smithRNA_fasta
        else mkdir $smithRNA_fasta
fi

if [ -f $main_outputs ]
        then rm $main_outputs
        else mkdir $main_outputs
fi

infoseq -only -length $home/$organism"_mit.fasta" > $home/length
len=$(cat $home/length | grep -v "Length")

#if shift option is activated (>1) the mito genome is cut at the selected position
if [[ $Shift -gt 1 && $Shift -le $len ]] ; then
        echo "SHIFT option activated. Mitochondrial genome cut at $Shift position" > $home/smith.log
        echo ""
        cat $home/smith.log 
	sleep 4
        seqret -send $(($Shift-1)) $home/$organism"_mit.fasta" -out $home/start.tmp.fa
        seqret -sbegin $Shift $home/$organism"_mit.fasta" -out $home/end.tmp.fa
        mv $home/$organism"_mit.fasta" $genome_mit/$organism"_old_mit.fasta"
        grep -v ">" $home/start.tmp.fa > $home/start.tmp
        cat $home/end.tmp.fa $home/start.tmp | sed '/^$/d' > $home/$organism"_mit.tmp.fasta"
        seqret $home/$organism"_mit.tmp.fasta" -outseq $home/$organism"_mit.fasta"
        rm $home/*tmp* 

fi

rm $home/length

#cut the first 1000 bps of the mito genome and paste at the end of the sequence (check warning)
seqret -send 1000 $home/$organism"_mit.fasta" -out $home/start.tmp.fa
seqret -sbegin 1001 $home/$organism"_mit.fasta" -out $home/end.tmp.fa

grep -v ">" $home/start.tmp.fa > $home/start.tmp

cat $home/end.tmp.fa $home/start.tmp | sed '/^$/d' > $home/$organism"_mit_cut.tmp.fasta"

seqret $home/$organism"_mit_cut.tmp.fasta" -outseq $home/$organism"_mit_cut.fasta"

rm $home/*tmp*

echo "" >> $home/smith.log
echo "START BUILDING BOWTIE2 INDEXES AT $(date +%X)" >> $home/smith.log
cat $home/smith.log | tail -n 2
sleep 10

#create bowtie indexes of the two genomes provided
bowtie2-build -f $home/$organism"_mit.fasta" $genome_mit/$organism"_mit" --threads $threads
bowtie2-build -f $home/$organism"_mit_cut.fasta" $genome_mit/$organism"_mit_cut" --threads $threads

cp $home/$organism"_mit.fasta" $genome_mit
mv $home/$organism"_mit_cut.fasta" $genome_mit

if [ -f $home/$organism"_nuc.fasta" ]; then
	bowtie2-build -f $home/$organism"_nuc.fasta" $genome_nuc/$organism"_nuc" --threads $threads
	cp $home/$organism"_nuc.fasta" $genome_nuc
fi


cd $fastq_folder;

##################################################################### TRIMMING ######################################################################################
#SE -> cutadapt
#PE -> fastp
#For paired-end (PE) trimming, the fastq files must follow this format: sampleX_1.fastq.gz, sampleX_2.fastq.gz. For single-end (SE) trimming, use the format

#trimming variables 
samples=$(ls *fastq.gz | sed "s/\_[1|2].fastq.gz//g" | sort | uniq);
num_threads_trim=5;


for sample in $samples; do 
	echo $sample >> $home/$organism"_samples.txt";	
done

sed -i '/^$/d' $home/$organism"_samples.txt"

echo "" >> $home/smith.log
echo ""

sleep 4

if [[ "$Trimming" == "NO" ]] ; then
	for sample in $samples; do
		cp $fastq_folder/$sample"_1.fastq.gz" $trimming
		mv $trimming/$sample"_1.fastq.gz" $trimming/$sample"_1.trimmed.fastq.gz"
		gunzip $trimming/$sample"_1.trimmed.fastq.gz"
	done

elif [[ "$Trimming" == "PE"  && "$adapter_sequence" != "NA" && "$adapter_sequenceR2" != "NA" ]] ; then

	echo "START TRIMMING AT $(date +%X)" >> $home/smith.log
	echo "Trimming performed in PE mode with adapters:" >> $home/smith.log
	echo "R1 adapter=$adapter_sequence" >> $home/smith.log
	echo "R2 adapter=$adapter_sequenceR2" >>$home/smith.log
	echo ""
	cat $home/smith.log | tail -n 4
	sleep 4

        for sample in $samples; do
                echo "fastp -c -x -g -a $adapter_sequence --adapter_sequence_r2 $adapter_sequenceR2 -w $num_threads_trim -i $fastq_folder/$sample"_1.fastq.gz" -I $fastq_folder/$sample"_2.fastq.gz" -o $trimming/$sample"_1.trimmed.fp.fastq" -O $trimming/$sample"_2.trimmed.fp.fastq"" >> $trimming/parallel.sh
        done

elif [[ "$Trimming" == "PE"  && "$adapter_sequence" != "NA"  && "$adapter_sequenceR2" == "NA" ]] ; then
	echo "START TRIMMING AT $(date +%X)" >> $home/smith.log
        echo "Trimming performed in PE mode with adapters:" >> $home/smith.log
        echo "R1 adapter=$adapter_sequence" >> $home/smith.log
        echo "R2 adapter=$adapter_sequenceR2" >>$home/smith.log
        echo ""
        cat $home/smith.log | tail -n 4
        sleep 4

        for sample in $samples; do
		echo "fastp -c -x -g -a $adapter_sequence -w $num_threads_trim -i $fastq_folder/$sample"_1.fastq.gz" -I $fastq_folder/$sample"_2.fastq.gz" -o $trimming/$sample"_1.trimmed.fp.fastq" -O $trimming/$sample"_2.trimmed.fp.fastq"" >> $trimming/parallel.sh
        done

elif [[ "$Trimming" == "PE"  && "$adapter_sequence" == "NA"  && "$adapter_sequenceR2" != "NA" ]] ; then
	echo "START TRIMMING AT $(date +%X)" >> $home/smith.log
        echo "Trimming performed in PE mode with adapters:" >> $home/smith.log
        echo "R1 adapter=$adapter_sequence" >> $home/smith.log
        echo "R2 adapter=$adapter_sequenceR2" >>$home/smith.log
        echo ""
        cat $home/smith.log | tail -n 4
        sleep 4

        for sample in $samples; do
                echo "fastp -c -x -g --adapter_sequence_r2 $adapter_sequenceR2 -w $num_threads_trim -i $fastq_folder/$sample"_1.fastq.gz" -I $fastq_folder/$sample"_2.fastq.gz" -o $trimming/$sample"_1.trimmed.fp.fastq" -O $trimming/$sample"_2.trimmed.fp.fastq"" >> $trimming/parallel.sh
        done


elif [[ "$Trimming" == "PE"  && "$adapter_sequence" == "NA"  && "$adapter_sequenceR2" == "NA" ]] ; then
	echo "START TRIMMING AT $(date +%X)" >> $home/smith.log	
	echo "Trimming performed in PE mode with adapter audotection" >> $home/smith.log
	echo "" >> $home/smith.log
	cat $home/smith.log | tail -n 3
	sleep 4

        for sample in $samples; do
		echo "fastp -c -x -g -w $num_threads_trim -i $fastq_folder/$sample"_1.fastq.gz" -I $fastq_folder/$sample"_2.fastq.gz" -o $trimming/$sample"_1.trimmed.fp.fastq" -O $trimming/$sample"_2.trimmed.fp.fastq"" >> $trimming/parallel.sh
        done

elif [[ "$Trimming" == "SE" && "$adapter_sequence" != "NA" ]] ; then
	echo "START TRIMMING AT $(date +%X)" >> $home/smith.log
	echo "Trimming performed in SE mode with adapter:" >> $home/smith.log
	echo "R1 adapter=$adapter_sequence" >> $home/smith.log
	echo "" >> $home/smith.log
	cat $home/smith.log | tail -n 4
	sleep 4

        for sample in $samples; do
                echo "cutadapt  --match-read-wildcards  --times 1 -e 0.1 -O 5 --cores $num_threads_trim --quality-cutoff 6 -m 15 --discard-untrimmed -a $adapter_sequence -o $trimming/$sample"_1.trimmed.ctdp.fastq" $fastq_folder/$sample"_1.fastq.gz"" >> $trimming/parallel.sh
       done
	
fi

cd $trimming

if [[ -f parallel.sh && $num_samples -gt 4 ]]; then
	parallel -j 4 < parallel.sh;
elif [[ -f parallel.sh && $num_samples -le 4 ]]; then
	parallel -j $num_samples < parallel.sh;
fi

rm -f *2.trimmed.fp.fastq


#Print the name of each sample within the header of the reads, after @ symbol
for trimmed in $(ls *.trimmed*); do
        a=$(echo $trimmed | sed s"/.trimmed.*//g");
        sed -i '1~4 s/^@/@'$a'_/g' $trimmed;
        echo "fastqc -t $num_threads_trim -o $fastqc $trimmed" >> parallel_fastqc.sh
	echo "gzip $trimmed" >> parallel_gzip.sh
done

#Quality check (fastqc)
if [ $num_samples -gt 4 ]; then
	parallel -j 4 < parallel_fastqc.sh;
else
	parallel -j $num_samples < parallel_fastqc.sh;
fi

#zip 
if [ $num_samples -gt 5 ]; then
        parallel -j 5 < parallel_gzip.sh;
else
        parallel -j $num_samples < parallel_gzip.sh;
fi

###################################################################ALIGNMENTS#################################################################
echo "" >> $home/smith.log
echo "START ALIGNMENTS AT $(date +%X)" >> $home/smith.log
echo ""
cat $home/smith.log | tail -n 2

for infile in $(ls *.trimmed* | sed "s/.trimmed.*//g"); do
        outfile=$infile\_on_MITO1.sam
        file=$infile\_MITO1_mapping.fastq
        outfile1=$infile\_on_MITO1.bam
        prefile1=$infile\_on_NUCL.sam
	prefile2=$infile\_on_NUCL.bam
        postfile1=$infile\_on_NUCL_unmapping.bam
        postfile2=$infile\_NUCL_unmapping.fastq
        postfile3=$infile\_on_MITO2.sam
        mito=$infile\_MitoUnique.bam


        #MITO1 Alignment: Alignment of total reads on mitichondrial genome
        bowtie2 -x $genome_mit/$organism"_mit" -U $infile*trimmed* -S $alignments/$outfile -N 1 -i C,1 -L 18 -p $threads --no-unal
        samtools sort -O BAM -o $alignments/$outfile1 $alignments/$outfile
        samtools flagstat $alignments/$outfile1 > $alignments/$infile.stats_MITO1.txt
        samtools depth -a $alignments/$outfile1 > $BedFiles/$outfile1.cov

	if [ -f $genome_nuc/$organism"_nuc.fasta" ]; then
		#NUCLEAR Alignment: Alignment of mito-mapping reads on nuclear genome
		samtools fastq  $alignments/$outfile1 > $alignments/$file
		bowtie2 -x $genome_nuc/$organism"_nuc" -q $alignments/$file -S $alignments/$prefile1 -i C,1 -L 22 -p $threads
		samtools sort -O BAM -o $alignments/$prefile2 $alignments/$prefile1
		samtools flagstat $alignments/$prefile2 > $alignments/$infile.stats_NUCL.txt
		
		#MITO2 Alignment: Alignment of mito-mapping and nuclear-non-mapping reads (mito-unique reads) on mitichondrial genome
		samtools view -b -f 4 $alignments/$prefile2 > $alignments/$postfile1
		samtools fastq $alignments/$postfile1 > $alignments/$postfile2
		bowtie2 -x $genome_mit/$organism"_mit" -q $alignments/$postfile2 -S $alignments/$postfile3 -N 1 -i C,1 -L 18 -p $threads
		samtools sort -O BAM -o  $alignments/$mito $alignments/$postfile3
		samtools flagstat $alignments/$mito > $alignments/$infile.stats_MITO2.txt
		samtools depth  -a $alignments/$mito > $BedFiles/$mito.cov

	else 
		mv $alignments/$outfile1 $alignments/$mito;
	fi


        #remove intermediate files
	rm -f $alignments/*sam $alignments/*MITO1.bam $alignments/*NUCL.bam $alignments/*unmapping.bam $alignments/*fastq
done;

cd $alignments;

#Create stats of MITO1, NUCLEAR and MITO2 alignments  
for i in $(ls *stats_MITO1.txt | sed 's/.stats_MITO1.txt//g'); do
        name=$i;
        mapped=$(grep "mapped" $i.stats_MITO1.txt | head -n 1);
        echo $name $mapped >> stats_MITO1.txt;
done;

if [ -f $genome_nuc/$organism"_nuc.fasta" ]; then
	for i in $(ls *stats_NUCL.txt | sed 's/.stats_NUCL.txt//g'); do
        	name=$i;
        	mapped=$(grep "mapped" $i.stats_NUCL.txt | head -n 1);
        	echo $name $mapped >> stats_NUCL.txt
	done;

	for i in $(ls *stats_MITO2.txt | sed 's/.stats_MITO2.txt//g'); do
	        name=$i;
        	mapped=$(grep "mapped" $i.stats_MITO2.txt | head -n 1);
        	echo $name $mapped >> stats_MITO2.txt;
	done;
fi

#Creation of bed files of MITO2 alignments
for infile in $(ls *Unique.bam | sed 's/.bam//g'); do
	bedtools bamtobed -i $infile.bam > $BedFiles/$infile.bed
	bedtools sort -i $BedFiles/$infile.bed > $BedFiles/$infile'_sort.bed'
        rm $BedFiles/*Unique.bed
done

cd $Clustering

###############################################################CLUSTERING###################################################################################
output=$organism"_results"
output_centroids=$output.centroids.fa
output_selected_centroids=$output.centroids.selected.fa
clusters_bed="5.2_"$output.clusters.bedfiles
clusters_fasta="5.1_"$output.clusters.fasta


#Clustering output folders creation
if [ -f $output".fa" ]
	then rm $output".fa"
fi

if [ -f $output_selected_centroids ]
	then rm $output_selected_centroids
fi

if [ -d $clusters_bed ]
	then rm -rf $clusters_bed
fi

mkdir $clusters_bed

if [ -d $clusters_fasta ]
	then rm -rf $clusters_fasta
fi

mkdir $clusters_fasta


#extraction of mito-unique reads of all the samples in a fasta file
for bam in $alignments/*Unique.bam; do 
        samtools fasta -n $bam >> $output".fa"
done

echo "" >> $home/smith.log
echo "START CLUSTERING (%ID=$clusteringID) AT $(date +%X)" >> $home/smith.log
echo ""
cat $home/smith.log | tail -n 2

#order sequences by abundance and print a sorted fasta file
seq_ordered=$(grep -v ">" $output".fa" | sort | uniq -c | sort -k 1 -r -h  | sed  's/^[ \t]*//' | sed 's/ /\t/g'  | cut -f 2)
for seq in $seq_ordered; do
        grep -B 1 -w "$seq" $output".fa" >> $output.sorted.fa;
done

rm $output".fa"

#formatting the fasta file
sed -i 's/--//g' $output.sorted.fa
sed -i '/^$/d'  $output.sorted.fa

#clustering 
vsearch --cluster_smallmem $output".sorted.fa" --usersort --threads 12 --centroids $output"_centroids.fa" --sizeout --clusterout_id --clusterout_sort --consout $output"_consensus.fa" --id $clusteringID --profile clusters_genomecov.txt --msaout alignment-consensu.txt --uc vsearch_table.txt --relabel_keep --minseqlength 15 --clusters clusterid;

######################################################################## CLUSTERING THRESHOLDS #################################################################################################
echo "" >> $home/smith.log
echo "START CLUSTERS FILTERING (STRINGENCY=$stringency) AT $(date +%X)"  >> $home/smith.log
echo ""
cat $home/smith.log | tail -n 2

mv clusterid* $clusters_fasta;

cd $clusters_fasta;

#number of samples replication
num_samples=$(cat $home/$organism"_samples.txt" | wc -l)

#creation of INFO.txt file with reports for each cluster the following information: clusterID, total reads number, number of reads belonging to each sample 
echo -e "CLUSTER\tCLUSTER_SIZE\t$(paste -sd "\t" $home/$organism"_samples.txt")" > INFO.txt

#create temporary files with total and per sample size information per cluster
for cluster in clusterid*; do
	cluster_size=$(grep ">" $cluster | wc -l)
        echo -e "$cluster\t$cluster_size" >> File.tmp;
        n=1
        for sample in $(cat $home/$organism"_samples.txt"); do
                grep ">" $cluster | grep  "$sample" | wc -l >> $n.$sample.file.tmp;
                n=$(($n+1))
done
done

#put the infos of temporary files in file INFO.txt
paste File.tmp  *file.tmp  >> INFO.txt


#Threshold1 (T1) definition: number of unique size clusters calculation and extraction of Xth cluster. Xth is selected accprding to stringency (S) parameter -> Xth = N clusters of uniqe size * (1 – S)
cut -f 2 INFO.txt | grep -v "CLUSTER" | sort -n -r| uniq > clusters_uniq.tmp
num_clusters_uniq=$(cut -f 2 INFO.txt | grep -v "CLUSTER" | sort | uniq | wc -l)
raw=$(echo "$num_clusters_uniq*(1-$stringency)" | bc | sed  "s/\..*//g")
sort -n -r clusters_uniq.tmp | head -n $raw | tail -n 1 > T1
echo ""
echo "Threshold1=$(cat T1)"

rm clusters_uniq.tmp


#Threshold2 (T2) definition: calculated as T1, but on the reads belonging to the single replicates 
for sample in $(cat $home/$organism"_samples.txt"); do
	#find the number corrisponding to each replicate within the table INFO.txt
        column=$(head -n 1 INFO.txt | sed -e 's/\t/\n/g' | grep -w -n "$sample" | sed s'/:/\t/g' | cut -f 1)
	cut -f $column INFO.txt | grep -v "$sample" | sort -n -r| uniq > $sample.clusters_uniq.tmp
        num_clusters_uniq_size=$(cut -f $column INFO.txt | grep -v "$sample" | sort | uniq | wc -l)
        raw=$(echo "$num_clusters_uniq_size*(1-$stringency)" | bc | sed  "s/\..*//g")
        sort -n -r $sample.clusters_uniq.tmp | head -n $raw | tail -n 1 > $sample.T2
        echo "Threshold2 for $sample=$(cat $sample.T2)"
done

rm *tmp

##########Check Reads mapping between the start and the end of the mito genome##############

#align all the reads on the cut  mito genome
bowtie2 -x $genome_mit/$organism"_mit_cut" -U $trimming/*trimmed*.fastq.gz -S $alignments/$organism"_mit_cut.sam" -N 1 -i C,1 -L 18 -p 10 --no-unal;

samtools sort -O BAM -o $alignments/$organism"_mit_cut.bam" $alignments/$organism"_mit_cut.sam";

samtools depth -a $alignments/$organism"_mit_cut.bam" > $BedFiles/$organism"_mit_cut.cov";


#max coverage calculation of the 10 upstairs and downstairs positions from the cutting point 
samtools faidx $genome_mit/$organism"_mit.fasta"
length=$(cut -f 2 $genome_mit/$organism"_mit.fasta.fai"); 
pos=$(($length - 1000)); 

max1=$(grep -B 10 -w "$pos" $BedFiles/$organism"_mit_cut.cov" | cut -f 3 | sort -h | tail -n 1)
max2=$(grep -A 10 -w "$pos" $BedFiles/$organism"_mit_cut.cov" | cut -f 3 | sort -h | tail -n 1)
T=$(cat T1)

#If max values exceed T1 then print the warning message 
if [[ $max1  -gt $T && $max2 -gt $T ]] ; then 
	echo "" >> $home/smith.log
	echo "Warning: some reads mapping through the end and the beginning of the mitochondrial genome are found. We suggest to use the shift option" >>  $home/smith.log; 
	cat $home/smith.log | tail -n 2;

fi

rm $alignments/$organism"_mit_cut.sam" $alignments/$organism"_mit_cut.bam" $BedFiles/$organism"_mit_cut.cov" $genome_mit/$organism"_mit_cut."*

#Print clusters passing or not T1 (FILTER_T1.txt) and T2 (FILTER_T2.txt)
for cluster in clusterid*; do
        Cluster_size=$(grep -c ">" $cluster)
        T1=$(cat T1)
        if [ $Cluster_size -ge $T1 ]
        then
        echo -e "$cluster\t$Cluster_size\tPASS" >> FILTER_T1.txt
        else
        echo -e "$cluster\t$Cluster_size\tNOT_PASS" >> FILTER_T1.txt
        fi
        for sample in $(cat $home/$organism"_samples.txt"); do
                Sample_size=$(grep $sample $cluster | wc -l )
                T2=$(cat $sample.T2)
                if [ $Sample_size -ge $T2 ]
                then
                        echo -e "$cluster\t$sample\t$T2\t$Sample_size\tPASS" >> replicates_threshold2.txt
                else
                        echo -e "$cluster\t$sample\t$T2\t$Sample_size\tNOT_PASS" >> replicates_threshold2.txt
                fi
        done

        #Count the number of replicates passing or not passing T2 
        N_rep_PASS=$(grep -w $cluster replicates_threshold2.txt | grep -w "PASS" | wc -l)
        N_rep_NOT_PASS=$(grep -w $cluster replicates_threshold2.txt | grep -w "NOT_PASS" | wc -l)

        #If the number of replicate passing T2 is >= M the cluster is PASS
        echo -e "$cluster\t$N_rep_PASS\t$N_rep_NOT_PASS" | awk -v num=$min_rep '{ if ($2>=num) print $1"\t"$2"\t"$3"\t""PASS"; else print $1"\t"$2"\t"$3"\t""NOT_PASS"}' >> FILTER_T2.txt

done

rm replicates_threshold2.txt

#formatting files
sed -i '1s/^/CLUSTER\tSIZE\tFILTER\n/' FILTER_T1.txt
sed -i '1s/^/CLUSTER\tN_REP_PASS\tN_REP_NOT_PASS\tFILTER\n/' FILTER_T2.txt

#print T1 results
echo ""
echo "CLUSTERING THRESHOLD 1:"

sleep 2

cat FILTER_T1.txt

#print T1 results
echo ""
echo "CLUSTERING THRESHOLD 2:"

sleep 2

cat FILTER_T2.txt


#Elimination of clasters NOT passing T1 e T2
for clusters_NOT_PASS in $(cat FILTER_T1.txt FILTER_T2.txt | grep -w "NOT_PASS" | cut -f 1 | sort | uniq); do
	rm $clusters_NOT_PASS;
done

#Renaming PASS clusters
for clusters_PASS in clusterid*; do
	mv $clusters_PASS ${clusters_PASS#clusterid}.fa; 
done


mv *T2 ../
mv T1 ../
mv FILTER* ../
mv INFO.txt ../

#Creation of BED files of "PASS" clusters
cp $BedFiles/*sort.bed .;
cat *bed > ALL.bed;
for cluster in $(ls *fa); do
        for seq in $(grep ">" $cluster | sed 's/>//g' | sed 's/;/\t/g' | cut -f 1); do 
		grep -w $seq ALL.bed >> $Clustering/$clusters_bed/${cluster%"fa"}"bed" ;
	done;
done

rm *bed

# *.genome file creation with name and length of the mithocondrial genome sequence. It will be used in the next step
name=$(cut -f 1 $genome_mit/$organism"_mit.fasta.fai")
length=$(cut -f 2 $genome_mit/$organism"_mit.fasta.fai")
echo -e $name"\t"$length > $Clustering/$organism".genome"


cd $Clustering/$clusters_bed

#Production of TOT, 3 and 5 coverage files (genomecov) for each cluster
for bed in *bed
        do
        bedtools sort -i $bed > ${bed%.bed}.sorted.bed
        bedtools genomecov -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome > ${bed%.bed}.genomecov
        bedtools genomecov -5 -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome > ${bed%.bed}.genomecov.5
        bedtools genomecov -3 -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome  > ${bed%.bed}.genomecov.3
done

rm $Clustering/*genome
cd $Clustering/$clusters_fasta

#fasta files creation of centroids belonging to PASS clusters
for i in $(ls *fa | sed 's/.fa//g'); do
grep -w -A 1 $(echo clusterid=$i) $Clustering/$output"_centroids.fa" > $Clustering/$i"_centroid.fa";
done

cd $Clustering

#For each PASS centroid: define ID and position and print it in the header 
for i in $(ls *_centroid.fa | sed 's/_centroid.fa//g'); do
        centroid=$(grep ">" $i"_centroid.fa" | sed 's/>//g' | sed 's/;/\t/g' | cut -f 1)
        position=$(grep -w $centroid $Clustering/$clusters_bed/*sorted.bed | cut -f 2,3,6 | sed 's/\t/;/g')
        sed -i "1 s/$/;$position/g" $i"_centroid.fa"
done

#concatenate the PASS centroids in a single multifasta file
cat *_centroid.fa > $output"_centroids.PASS.fasta"

#multifasta formatting
sed -i 's/;/_/g' $output"_centroids.PASS.fasta" 
sed -i  's/>.*clusterid/>clusterid/g' $output"_centroids.PASS.fasta"
sed -i 's/=//g' $output"_centroids.PASS.fasta"
sed -i 's/_/_pos/2' $output"_centroids.PASS.fasta"
sed -i 's/_/_strand/4' $output"_centroids.PASS.fasta"

rm *_centroid.fa

######################################################### PLOTS ##################################################
echo "" >> $home/smith.log
echo "START PLOTTING DATA AT $(date +%X)"  >> $home/smith.log
echo ""
cat $home/smith.log | tail -n 2

cd $Plots

#copy files for plotting
cp $Clustering/T1 .
cp $Clustering/*T2 .
cp $Clustering/$output"_centroids.PASS.fasta" .


#copy and formatting coverage files in tables COV1 and COV2 reads coverage of MITO1 and MITO2 alignments
if [ -f $genome_nuc/$organism"_nuc.fasta" ]; then
	cp $BedFiles/*Unique.bam.cov .
fi

cp $BedFiles/*MITO1.bam.cov .

#Position file creation
for i in *MITO1.bam.cov; do
        cut -f 2 $i > POS.temp;
done

#Creatoion of temporary files of MAP1 coverage for each sample
for i in *MITO1.bam.cov; do
        cut -f 3 $i > $i.COV1.temp;
done

#Merge COV1 files and calculation of total coverage for each position
paste *COV1.temp | awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'> SUM1.temp

#COV1 table creation: coverage of TOTAL reads on mtochondrial genome
paste POS.temp *COV1.temp SUM1.temp > COV1.txt

if [ -f $genome_nuc/$organism"_nuc.fasta" ]; then

        #Creatoion of temporary files of MAP2 coverage for each sample
        for i in *MitoUnique.bam.cov; do
                cut -f 3 $i > $i.COV2.temp;
        done

        #Merge COV2 files and calculation of total coverage for each position
        paste *COV2.temp | awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'> SUM2.temp

        #COV2 table creatio: coverage of mitochondrial unique reads on mitochondrial genome
        paste POS.temp *COV2.temp SUM2.temp > COV2.txt

        #Rscript execution coverage of total and mito-unique reads on mithocondrial genome
        Rscript $script1 COV1.txt COV2.txt $PWD $output"_centroids.PASS.fasta"
else
	#If nuclear genome is not provided we simulate a COV2 file which is identical to the COV 1 file 
	cp COV1.txt COV2.txt

	Rscript $script1 COV1.txt COV2.txt $PWD $output"_centroids.PASS.fasta"
fi

#Python script execution -> Coverage of each PASS cluster
$script2 -d $Clustering/$clusters_bed

rm *temp
rm *cov
rm T1 *T2 *_centroids.PASS.fasta

cd $home

################SMITHS FASTA FILE GENERATION#################################
echo "" >> $home/smith.log
echo "SMITHS FASTA FILE GENERATION AT $(date +%X)"  >> $home/smith.log
echo ""
cat $home/smith.log | tail -n 2


cp $Clustering/$output"_centroids.PASS.fasta" $smithRNA_fasta


#check if smiths in reverse are present and do the reverse complement
if [ $(grep ">" $smithRNA_fasta/$output"_centroids.PASS.fasta" | grep "strand-" | wc -l) -gt 0 ];

then
        cd $smithRNA_fasta

        #from multifa to single fasta
        seqret -ossingle -sequence $output"_centroids.PASS.fasta" -outseq seqoutall

        #list of reverse smith
        ls * | grep "strand-" > list_reverse

        sed -i 's/.fasta//g' list_reverse
        sed -i 's/strand-//g' list_reverse

        #for each smith in reverse do the reverse complement
	#In the header the strand is kept minus to maintain the information of the original orientation 
        for i in $(cat list_reverse); do
        revseq $i"strand-.fasta" -tag FALSE -outseq $i"strand+.fasta";
        done

        cat *+.fasta > presumptive_smithRNAs.fa

        rm *fasta list_reverse

else
        mv $smithRNA_fasta/$output"_centroids.PASS.fasta" $smithRNA_fasta/presumptive_smithRNAs.fa


fi

cd $home

#copy main outputs in the folder
mkdir $main_outputs/Plots
cp  $smithRNA_fasta/presumptive_smithRNAs.fa $main_outputs
cp $Plots/COV2.pdf $main_outputs/Plots
cp $Plots/COV2.clusters.pdf $main_outputs/Plots
cp $Plots/Plot.html $main_outputs/Plots
cp $Plots/COV2.replicates.pdf $main_outputs/Plots
cp $Plots/COV2.stats $main_outputs/Plots

#conda enviroment deactivation
conda deactivate

exit 0
