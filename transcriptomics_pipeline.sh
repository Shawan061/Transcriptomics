#!/bin/bash

#SBATCH --job-name=transcriptomics_pipeline_test      ## Job name for identification 
#SBATCH --partition=cpu-short     ## Partition (queue) to run in (e.g., cpu-short, gpu)       
#SBATCH --nodes=1                              ## Number of nodes to allocate 
#SBATCH --ntasks=1                             ## Number of tasks (often 1 for non-MPI jobs) 
#SBATCH --cpus-per-task=31            ## Number of CPU cores per task 
#SBATCH --mem=50G                    ## Memory allocation per node (e.g., 32 GB) 
#SBATCH --time=0-24:00:00             ## Time limit (D-HH:MM:SS)
#SBATCH --output=transcriptomics_pipeline_shawan_%j.log ## Standard output and error log file (%j = Job ID) 
 
 
# --- Print job information --- 
echo "Job started on $(date)" 
echo "Job ID: $SLURM_JOB_ID" 
echo "Running on nodes: $SLURM_JOB_NODELIST" 
echo "Working directory: $(pwd)" 


# --- Loading required modules ---
module load fastqc-0.12.0
module load trimmomatic-0.39
module load star-2.7.11
module load trimgalore-0.6.10
module load subread-2.1.1


cd /home/d_parashar/shawan/dengue_transcriptomics_project/fastq

#organizing the fastq files into folders
echo "                              <<<<<<<<<<<<<<<<<<<<<<<<<<<< making folders and moving >>>>>>>>>>>>>>>>>>>>>>>>>>>>"
file=$(ls | sed 's/_[12]\.fastq\.gz$//')
for f in $file; do
    mkdir ${f}
    mv ${f}_1.fastq.gz ${f} >/dev/null 2>&1
    mv ${f}_2.fastq.gz ${f} >/dev/null 2>&1
    echo "                **************** file moving done, going to next folder****************"
done

echo -e "\nFolders created and files moved successfully 👍\n"


echo -e "\n                              <<<<<<<<<<<<<<<<<<<<<<<<<<<< running fastqc >>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo -e "\ncurrently in: $main\n"
for f in $folders; do 
    echo -e "\n                          **************** Entering folder: $f ****************\n"
    cd $f || { echo "Failed to enter $f"; exit 1; }
    pwd
    mkdir -p fastqc_output
    echo -e "\n                ######################## running fastqc ########################\n"
    fastqc -t 8 *.fastq.gz -o fastqc_output/
    echo -e "\n                           **************** ✅fastqc for $f is done✅ ****************\n"
    cd "$main"
done

echo -e "\nfastqc is done for every file 👍\n"

#trimming


echo -e "\n                <<<<<<<<<<<<<<<<<<<<<<<<<<<< Trimming adaptors using trimgalore >>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"

folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo -e "\ncurrently in: $main\n"
for f in $folders; do 
	cd $f/fastqc_output
    echo -e "################### Processing folder: ${f} ###################\n"
	#unzipping the fastq output
	for d in *fastqc.zip; do
        unzip -oq "$d" | echo -e "unzip done for $d"
    done
    adaptor_count=0
    inner_pwd=$(pwd)
	inner_folders=$(ls -d */ | sed 's:/$::')
	for F in $inner_folders; do
		cd $F
        echo -e "\ncurrently in: $F"
		adap=$(grep -E '>>Adapter Content' fastqc_data.txt | awk '{print $3}')
        
        if [ "$adap" = "pass" ]; then
            echo -e "\n adaptor contamination status: $adap ✅\n<<No adaptor detected in ${F}>> \n "
		else 
            echo -e "\n adaptor contamination status: $adap 🚫\n<<Adaptor detected in ${F}>> \n "
            adaptor_count=$((adaptor_count + 1))
        fi
        cd "$inner_pwd"
    done
    echo -e "Adaptor count in ${f}: $adaptor_count \n"
    cd ..
    if [ "$adaptor_count" -gt 0 ]; then
        echo -e "\nAdaptor detected in atleast one of the files in ${f}, proceeding to trimgalore ✂️\n"
        echo -e "\n                ######################## running trim_galore ########################\n"
        trim_galore --paired ${f}_1.fastq.gz ${f}_2.fastq.gz -j 8 -o trimgalore_output
        cd trimgalore_output || { echo "Failed to enter trimgalore_output"; exit 1; }
        mv ${f}_1_val_1.fq.gz ${f}_trimgalore_1.fastq.gz
        mv ${f}_2_val_2.fq.gz ${f}_trimgalore_2.fastq.gz
    else
        echo -e "No adaptor detected in any of the files in ${f}, Skipping trimgalore ⏭️\n"
        mkdir -p trimgalore_output
        cp *.fastq.gz trimgalore_output/
        cd trimgalore_output
        mv ${f}_1.fastq.gz ${f}_trimgalore_1.fastq.gz
        mv ${f}_2.fastq.gz ${f}_trimgalore_2.fastq.gz
    fi
    echo -e "\n done with ${f} trimming, moving to next folder ✅\n"
    cd "$main"
done

echo -e "\nAdaptor trimming is done for every file 👍\n"

#activaing conda environment for komplexity
source activate /home/d_parashar/.conda/envs/complexity_env




echo -e "\n                      <<<<<<<<<<<<<<<<<<<<<<<<<<<< making overrepresented file >>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo -e "\nall the folders are present in: $main\n"
for f in $folders; do
	echo -e "\n                          **************** processing on sample: $f ****************\n"
	cd $f || { echo "Failed to enter $f"; exit 1; }
	pwd
	cd fastqc_output
	pwd
	inner_folders=$(ls -d */ | sed 's:/$::')
	for F in $inner_folders; do
		cd ${F}
		#filter out all the overrepresented sequences
		grep -E '^[ATGCN]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]*' fastqc_data.txt | awk '{print ">"$4"_"NR "\n"$1 }' > result.fasta
		#running komplexity to identify low complex reads
		kz --fasta < result.fasta > out.txt
		#filter which reads are no hit
		grep -A1 '^>[^No]' result.fasta | grep -v '^--$' > primer.fasta
		#filter out low complexity reads
		awk '$4 < 0.4 { print $1 }' out.txt | while read id; do grep -wA1 ">$id" result.fasta; done > low_complex.fasta
		cat low_complex.fasta primer.fasta > trim.fasta
        cp trim.fasta "$main/$f/trimgalore_output/"
		cd ..
	done	
	cd "$main"
	echo -e "\nall necessary contamination reads for sample $f are saved in trim.fasta file"	
done

#deactivating the conda environment
conda deactivate

echo -e "\nTrimable overrepresented sequence has been made for every file 👍\n"




echo "           <<<<<<<<<<<<<<<<<<<<<<<<<<<< Removing overrepresented sequences using trimmomatic >>>>>>>>>>>>>>>>>>>>>>>>>>>>"
folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo "all the folders are present in: $main"
for f in $folders; do 
    cd $f/trimgalore_output || { echo "Failed to enter $f/trimgalore_output"; exit 1; }
    echo -e "\n                          **************** Entering folder: $f ****************\n"
    pwd
    mkdir -p trimmomatic_output
    echo -e "\n                ######################## running Trimmomatic ########################\n"
    if [ ! -s  trim.fasta ]; then
        echo -e "\ntrim.fasta is empty, skipping trimmomatic for ${f} ⏭️\n"
        mv ${f}_trimgalore_1.fastq.gz ${f}_trimmomatic_1.fastq.gz
        mv ${f}_trimgalore_2.fastq.gz ${f}_trimmomatic_2.fastq.gz
        cp ${f}_trimmomatic_1.fastq.gz ${f}_trimmomatic_2.fastq.gz trimmomatic_output/
        cd trimmomatic_output
        mkdir -p fastqc_output
        fastqc -t 10 *.fastq.gz -o fastqc_output/
    else
        echo -e "\ntrim.fasta file found, runing trimmomatic on ${f} ✂️\n"
        trimmomatic PE -threads 8 -phred33 ${f}_trimgalore_1.fastq.gz ${f}_trimgalore_2.fastq.gz ${f}_trimmomatic_paired_1.fastq.gz ${f}_trimmomatic_unpaired_1.fastq.gz ${f}_trimmomatic_paired_2.fastq.gz ${f}_trimmomatic_unpaired_2.fastq.gz ILLUMINACLIP:trim.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        mv "${f}_trimmomatic_paired_1.fastq.gz" "${f}_trimmomatic_paired_2.fastq.gz" trimmomatic_output
        cd trimmomatic_output
        mv ${f}_trimmomatic_paired_1.fastq.gz ${f}_trimmomatic_1.fastq.gz
        mv ${f}_trimmomatic_paired_2.fastq.gz ${f}_trimmomatic_2.fastq.gz
        mkdir -p fastqc_output
        fastqc -t 10 *.fastq.gz -o fastqc_output/
    fi
    cd "$main"
done

echo -e "\nTrimable overrepresented sequence has been removed from every file 👍\n"


ulimit -n 100000
#creating BAM files from the trimmed .fq files using STAR

echo -e "\n                              <<<<<<<<<<<<<<<<<<<<<<<<<<<< STAR Alignment >>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"
folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo "\nall the folders are present in: $main\n"
for f in $folders; do
    cd $f/trimgalore_output/trimmomatic_output || { echo "Failed to enter ${f} — skipping to next."; cd "$main"; continue; }
    echo -e "\n                           **************** running process on ${f} ****************\n"
    pwd
    mkdir -p "bam_files"
    echo -e "\n                           **************** running STAR ****************\n"
    STAR --runThreadN 30 --genomeDir /home/d_parashar/shawan/reference_genome --readFilesIn *_1.fastq.gz *_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${f}_ --outSAMtype BAM SortedByCoordinate
    echo -e "\n                           **************** BAM file for ${f} has been generated ✅ ****************\n"
    mv ${f}_Aligned.sortedByCoord.out.bam ${f}.bam
    mv ${f}.bam bam_files/
    cd "$main"
done
echo "\nall BAM files has been generated 👍\n"


#creating count matrix using featureCounts

echo "                               <<<<<<<<<<<<<<<<<<<<<<<<<<<< Featurecounts >>>>>>>>>>>>>>>>>>>>>>>>>>>>"
folders=$(ls -d */ | sed 's:/$::')
main=$(pwd)
echo "all the folders are present in: $main"
mkdir featureCounts
for f in $folders; do
    cd ${f}/trimgalore_output/trimmomatic_output/bam_files || { echo "Failed to enter $f/trimmed/trimmomatic_output — skipping to next."; cd "$main"; continue; }
    echo "processing on ${f}"
    cp ${f}.bam ${main}/featureCounts
    echo "successfully copied ${f} to ${main}/featureCounts"
    cd "$main"
done
echo -e "\nfile copy done, starting featureCounts\n"
cd ${main}/featureCounts
featureCounts -a /home/d_parashar/shawan/reference_genome/gencode.v48.chr_patch_hapl_scaff.annotation.gtf -o counts.txt -g gene_name -T 30 -p *.bam

echo -e "\nFeatureCounts has been run successfully, counts.txt file is generated in ${main}/featureCounts 👍\n"

#srun --partition=cpu-debug --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=100G --time=0-2:00:00 bash full_pipeline_copy.sh
