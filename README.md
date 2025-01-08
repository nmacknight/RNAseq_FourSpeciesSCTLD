# RNAseq Bioinformatics Pipeline

**Project:** ***Acropora cervicornis***, ***Porites astreoides***, ***Montastraea cavernosa***, and ***Orbicella faveolata*** were exposed to Stony Coral Tissue Loss Disease. Tissue samples were collected to investigate the gene expression of each component of the coral holobiont with particular consideration of the Bacteria Gene expression as this focus has never been investigated and may lead to results that identify a bacteria disease causing agent and mechanistic etiology. 

**Goal:** Process sequence reads to build Metatranscriptomes, isolate holobiont components (Host, Algal Symbiont, Bacteria), identify single copy Orthologs, and annotate those Orthologs for statistical analysis in R.  

> [!TIP]
> Click the arrows to expand each section

<details>

<summary>Downloading Software</summary>

# Downloading Software

> Notes for myself: I am installing this software on a server which required me to install basic software into my user folder because it was a blank slate. I did find a parent folder that had a lot of software, but there were inconsistencies on whether I could execute using those pre-installed software tools and certainly could not write new software into that parent folder. Also, some software, I believe Trinity and Samtools were successfully compiled either by ethernet connection in the aoml building like how I mentioned below, and/or by being in a conda environment and then performing the compiling or running the code.

> Notes for you, another user: I had user permission limitations due to this being conducted on a government server and did not and could not be granted admin privileges or the ability to use sudo. So, while this is how I installed and compiled the necessary software, if you have admin priviliges on your server, you are likely to run into less errors, hopefully. It also means the exact code to install the software may differ for you. 

### Anaconda
Download script to install anaconds:
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
```

Make the script executable:
```
chmod +x Anaconda3-2021.05-Linux-x86_64.sh
```

Execute Script (installs):
```
./Anaconda3-2021.05-Linux-x86_64.sh
```

Verify Installation:
```
anaconda --version
python --version
```

### FASTP
> Fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
```
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

### TRINITY
> Trinity is a de novo transcriptome assembler for RNA-Seq data.

Requires dependencies to be installed too. 
/home/cns.local/nicholas.macknight/software
Trinity Dependencies: Jellyfish, salmon, bowtie, bowtie2, samtools (which requires ncurses development tool dependency)

**Install jellyfish**
```
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-linux
chmod +x jellyfish-linux
mv jellyfish-linux jellyfish
```

**Install salmon**
```
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz
tar -xzvf salmon-1.5.2_linux_x86_64.tar.gz
mv salmon-1.5.2_linux_x86_64 salmon
```

**Install bowtie**
```
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip
unzip bowtie-1.3.0-linux-x86_64.zip
mv bowtie-1.3.0 bowtie
```

**Install bowtie2**
```
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip
unzip bowtie2-2.4.4-linux-x86_64.zip
mv bowtie2-2.4.4 bowtie2
```

**Install curses development files for Debian/Ubuntu**
```
wget http://ftp.us.debian.org/debian/pool/main/n/ncurses/ncurses_6.2+20201114.orig.tar.gz
tar -xzvf ncurses_6.2+20201114.orig.tar.gz
cd ncurses-6.2+20201114
./configure --prefix=/home/cns.local/nicholas.macknight/software
make
make install
cd ..
```

> Could not get samtools to identify ncurses during the samtools configure stage. So I ran ./configure in the samtools folder like: "./configure --without-curses" curses is required for the tview command in samtools. we will see if this is necessary later one.

**Install samtools**
```
wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2
tar -xjvf samtools-1.13.tar.bz2
cd samtools-1.13
./configure --without-curses --prefix=/home/cns.local/nicholas.macknight/software/samtools-1.13
make
make install
cd ..
```

**Install java (if not already installed)**
/home/cns.local/nicholas.macknight/software/
```
curl -s "https://get.sdkman.io" | bash
source "/home/cns.local/nicholas.macknight/.sdkman/bin/sdkman-init.sh"
sdk install java
```
**Need to install cmake for Trinity**
```
wget https://github.com/Kitware/CMake/releases/download/v3.28.3/cmake-3.28.3.tar.gz
tar -xzvf cmake-3.28.3.tar.gz
cd cmake-3.28.3
mkdir build

export PATH="/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/cmake-3.28.3/bin:$PATH"
source ~/.bashrc  # or source ~/.zshrc if you are using Zsh
```

**Install Trinity**
```
mkdir ~/Trinity
cd ~/Trinity
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.0/trinityrnaseq-v2.15.0.FULL.tar.gz
tar -xzvf trinityrnaseq-v2.15.0.FULL.tar.gz
cd trinityrnaseq-v2.15.0
make
```
> could not get make to work. #SUCCESS: ran the "make" command in the aoml building while connected to their ethernet internet. worked immediately. wasnt working at home due to user permission error. Fun.

Test if Trinity was installed properly
```
export TRINITY_HOME=/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0
cd sample_data/test_Trinity_Assembly/
./runMe.sh
```

**ParaFly**
```
wget https://github.com/ParaFly/ParaFly
```

Add Trinity and its dependencies to PATH
```
export TRINITY_HOME=/space/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0
export PATH=$TRINITY_HOME:$TRINITY_HOME/util:$TRINITY_HOME/Inchworm:$TRINITY_HOME/Chrysalis:$TRINITY_HOME/Butterfly:$PATH
```

### BLAST 
> Extract Coral Only (or Symbiont Only or Bacteria Only) Transcripts
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
```

### cd-hit
> Create BLAST Database
```
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
tar -xvzf cd-hit-v4.8.1-2019-0228.tar.gz
make
cd cd-hit-auxtools
make
```

### cdbfasta
> CDBFASTA is a tool used for indexing and querying large sequence sets in the FASTA format efficiently. It creates a hash index of sequences, allowing for fast retrieval and comparison of sequences based on their headers or identifiers.
```
git clone https://github.com/gpertea/cdbfasta.git
cd cdbfasta
make
```
### UNIPROT
> This is the database that will be used to annotate sequences and assign an Entry ID and an E-value number which is a score on how well the sequence matches the database sequence. Low E values represent high sequence matching, so you want low e-values. As a standard bench mark, E-values with less than e^-05 or 0.00001 are retained and anything greater than that is filtered out.
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```

### BOOST
> This may have not been utilized. 
```
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
```

### TOPHAT
```
git clone https://github.com/infphilo/tophat.git
./bootstrap
./configure
```

### DB_File
```
cpan
install DB_File
exit
```

### SAMTOOLS
```
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar -xjf samtools-1.14.tar.bz2
./configure --without-curses
make
```
> ERROR : Permission denied. Probably fixed by trying once connected to aoml ethernet


### BBMAP
```
wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz
tar -xvzf BBMap_38.90.tar.gz
```

### SALMON
```
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/Salmon-1.5.2_linux_x86_64.tar.gz
tar -xvzf Salmon-1.5.2_linux_x86_64.tar.gz
```

### OrthoFinder
```
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder_source.tar.gz
tar xzf OrthoFinder_source.tar.gz
python OrthoFinder_source/orthofinder.py -h # to verify installation
```

### Transdecoder
```
wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz
tar -xvzf TransDecoder-v5.7.1.tar.gz
```

### Singularity
```
wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/archived/trinityrnaseq.v2.15.0.simg
tar -xzf trinityrnaseq.v2.15.0.simg
cd singularity
./mconfig
make -C builddir
```

### Libevent
```
wget http://monkey.org/~provos/libevent-1.4.14b-stable.tar.gz
tar xzf libevent-1.4.14b-stable.tar.gz
cd libevent-1.4.14b-stable
./configure --prefix=/opt/libevent
make
make install #permission denied
```

### ncurses
```
wget http://invisible-mirror.net/archives/ncurses/current/ncurses.tar.gz
tar xzf ncurses.tar.gz
cd ncurses-6.4-20240224
./configure
make
```

### tmux
```
wget https://github.com/tmux/tmux/releases/download/3.4/tmux-3.4.tar.gz
tar xvf tmux-3.4.tar.gz
cd tmux-3.4/
./configure # needs libevent. Failed. 
make
```

### BLAST
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.15.0+-x64-linux.tar.gz
cd ncbi-blast-2.15.0+
export PATH=$PATH:/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin
```
### BUSCO
> BUSCO is a tool to determine how "complete" a transcriptome is by identiying the presence of highly conserved genes in your assembled transcriptome. It is a way to make sure even though you may have a lot of transcripts that it represents a "complete" biological organism and a nice checkpoint that things are on the right track.
> Now, BUSCO requires a lot of dependencies or user permissions that I could not get (at the time of writing this).
Below are the resources I have in my attempts to install BUSCO. 
```
TRINITY QUALITY ASSESSMENT
BUSCO

git clone https://gitlab.com/ezlab/busco.git
cd busco/
export PATH=:/home/cns.local/nicholas.macknight/software/busco/bin$PATH

# Install using conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b -p /home/cns.local/nicholas.macknight/software/busco/miniconda
eval "$(/home/cns.local/nicholas.macknight/software/busco/miniconda/bin/conda shell.bash hook)"
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda update --all

conda activate 
conda create -n busco_env
conda activate busco_env
conda install busco=4.0.6


conda install -c conda-forge mamba

# EXAMPLE of busco command:
busco -i [SEQUENCE_FILE] -m [MODE] [OTHER OPTIONS]

# Modified for our test data:
~/.local/bin/busco -i /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -m transcriptome
~/.local/bin/busco -i /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/final_coral_reference_transcriptome.fa -m transcriptome
```
</details>


<details>
<summary>Downloading Data to Server</summary>

# Downloading Data to Server


> [!TIP]
> Moving Sequence data from the sequencer to your server/local computer is going to be project specific. My first recommendation is to have a discussion with your sequencer on how to make the necessary transfer of data. Hopefully they have a tutorial with clear instructions. In the past I have used filezilla. For this project they uploaded the data on box but the amount of data was too large for box to be able to download to my computer or server. So, following my own advice I reached out to the sequencer and they transferred the data to an amazon host server and gave me instrucitons on how to transfer the data to my local computer, where I redundantly saved it to a local external hard drive and then moved it to my server. 

Took 12 hours to download all 1232 files. 
To Confirm all 1232 files were downloaded, use this command:
```
find ./sRosales_OfavSCTLD/ -maxdepth 1 -type f | wc -l
```
It will tell you how many unique files are in the folder sRosales_OfavSCTLD/
</details>

<details>
<summary>Fastp - PreProcessing</summary>
# Fastp - PreProcessing
> The purpose of Fastp is to remove adapters and low quality reads. 

export PATH=/home/cns.local/nicholas.macknight/software:$PATH; for file1 in ./RNARawData/sRosales_OfavSCTLD/*_R1_*.fastq.gz; do base=$(basename "$file1" _R1_001.fastq.gz); file2="./RNARawData/sRosales_OfavSCTLD/${base}_R2_001.fastq.gz"; fastp --verbose -i "$file1" -I "$file2" -o "./Fastp_ProcessedData/${base}_clean_R1.fastq.gz" -O "./Fastp_ProcessedData/${base}_clean_R2.fastq.gz"; done

# You can occassionally check on how many files are completed (after opening another terminal window) with:
```
find ./Fastp_ProcessedData/ -maxdepth 1 -type f | wc -l
```

###Concatenate Cleaned Fastp Files:

Created a script called **FastpMerge.sh** in /home/cns.local/nicholas.macknight/SCTLDRNA/ and is as follows:

```
# Input directory for fastq files
input_dir="/home/cns.local/nicholas.macknight/SCTLDRNA/Fastp_ProcessedData"

# Output directory for merged files
output_dir="/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/"

# Get a list of unique names by extracting the first two portions from the file names
names=($(ls "${input_dir}" | cut -d'_' -f1-2 | sort -u))

# Loop through each name in the list
for name in "${names[@]}"; do
    cat ${input_dir}/${name}_L00?_clean_R1.fastq.gz > "${output_dir}/${name}_R1_clean_merged.fastq.gz"
    cat ${input_dir}/${name}_L00?_clean_R2.fastq.gz > "${output_dir}/${name}_R2_clean_merged.fastq.gz"
done
```

Here are a few manual examples of simply concatenating the runs (Lane 1-8) of each sample: 

```
cat Ofav17-7Control-40_S18_L00?_clean_R1.fastq.gz > ../MergedFastpProcessedData/Ofav17-7Control-40_S18_R1_clean_merged.fastq.gz
cat Ofav17-7Disease-9_S23_L00?_clean_R1.fastq.gz > ../MergedFastpProcessedData/Ofav17-7Disease-9_S23_R1_clean_merged.fastq.gz
```
You can either modify a for lopp to do this or run it manually. 
</details>

<details>

<summary>Trinity</summary>
#Trinity
> Assembles Transcript sequences into de novo Transcriptomes

The most fundamental EXAMPLE of the code:
```
Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G
```

Modifying for our real data: (single lane of a single sample forward and reverse reads - this is a test that Trinity is run)
```
Trinity --seqType fq --left PastPA3Disease-5_S62_L001_clean_R1.fastq.gz --right PastPA3Disease-5_S62_L001_clean_R2.fastq.gz --CPU 64 --max_memory 400G 
```

Determine Memory availability of the server:
```
free -h #487Gb
```

Determine CPUs on the server:
```
lscpu #64
```


**Estimate Memory Usage per Sample:**

Trinity estimates memory usage based on the number of reads and k-mer content in the dataset. You can use the Trinity script util/insilico_read_normalization.pl to estimate memory requirements for each sample.

Example:
```
/path/to/trinity/util/insilico_read_normalization.pl --seqType fq \
  --JM 100G --left sample_R1.fastq --right sample_R2.fastq
  
/space/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/insilico_read_normalization.pl --seqType fq \
  --JM 100G --left PastPA3Disease-5_S62_L001_clean_R1.fastq.gz --right PastPA3Disease-5_S62_L001_clean_R2.fastq.gz 
```
Adjust the --JM parameter until it fits within your available memory.

Total Memory Requirement:
Once you have estimated the memory requirement for one sample, multiply it by the number of samples you want to process concurrently.
Total Memory = Memory per Sample * Number of Samples





all:
    	mkdir -p build
        cd build && cmake -DCMAKE_INSTALL_PREFIX="" ../ && make DESTDIR=../ install

clean:
      	@echo cleaning
        (cd build && make clean) || :
        rm -rf ./build ./bin




Test Trinity installation in /home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/sample_data/test_Trinity_Assembly with:
```
./runMe.sh
```
To successfully run Trinity, we need to tell the server where to look to find Trinity's dependencies. 
```
export TRINITY_HOME=/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0
export PATH=/home/cns.local/nicholas.macknight/software:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/samtools-1.14:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/bowtie:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/bowtie:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/bowtie2:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin:$PATH
export PATH=/home/cns.local/nicholas.macknight/software/anaconda3/bin:$PATH
export PATH=/home/cns.local/nicholas.macknight/anaconda3/:$PATH
export PATH=:/usr/bin/screen$PATH
export PATH=/bin:$PATH
export PATH=/bin/anaconda3/bin:$PATH
```
If you need to activate the conda environmnet:
```
source /bin/anaconda3/etc/profile.d/conda.sh
conda activate #activates conda environment
```
Here is an example of a software that needed to be installed once you have activated a conda environment:
```
conda install numpy # install numpy
cd ../../../home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/sample_data/test_Trinity_Assembly/
./runMe.sh # should work successfully
```

Final Trinity assemblies are written to /space/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/sample_data/test_Trinity_Assembly/trinity_out_dir.Trinity.fasta
[Info on interpreting Trinity Output](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Output-of-Trinity-Assembly)



count the number of transcripts
```
grep -c '^>' trinity_out_dir.Trinity.fasta
```

## Running Trinity 

### Porites astreoides

```
du -ch Past*R1* # 67gb
du -ch Past*R2* # 65gb
```

 Input data
  Left.fasta    30229 MByte
  Right.fasta   30229 MByte
  Number of unique KMERs: 3046076793
  Number of reads:        328328984 Output data
  Trinity.fasta 0 MByte

Runtime
=======
Start:       Sat Mar 23 12:41:39 EDT 2024
End:         Wed Mar 27 01:25:19 EDT 2024
Trinity   305020 seconds (84.73 hours or 3.5 days)
  Inchworm (phase 1 - read clustering)  22438 seconds
  Chrysalis (phase 1 - read clustering) 180607 seconds
  Rest (phase 2 - parallel assembly)	   101975 seconds
  
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/Trinity --seqType fq --left \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-5_S62_clean_R1_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-5-2_S72_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-15_S73_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Control-29_S75_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Pastpa3Control-24_S65_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Control-22_S74_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-8_S63_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-2_S61_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-12_S70_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-39_S77_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-34_S68_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-28_S67_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-4_S71_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-16_S69_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-10_S64_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Control-31_S76_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Control-21_S66_R1_clean_merged.fastq.gz \
--right /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-5_S62_clean_R2_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-5-2_S72_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Disease-15_S73_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Control-29_S75_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Pastpa3Control-24_S65_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA3Control-22_S74_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-8_S63_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-2_S61_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Disease-12_S70_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-39_S77_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-34_S68_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA2Control-28_S67_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-4_S71_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-16_S69_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Disease-10_S64_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Control-21_S66_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/PastPA1Control-31_S76_R2_clean_merged.fastq.gz \
 --CPU 64 --max_memory 400G --output trinity_out_dir_AllPastSamples_Lane1-8
 ```
 

### Acropora cervicornis 

Estimate total file size to estimate runtime.
```
du -ch Acer*R1* # 59G
du -ch Acer*R2* # 58G
```

```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/Trinity --seqType fq --left \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-31_S84_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-38_S89_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-12_S81_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-19_S83_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-7_S80_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-35_S85_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-39_S87_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-14_S90_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-4_S78_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-37_S86_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-40_S88_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-17_S82_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-6_S79_R1_clean_merged.fastq.gz \
--right \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-31_S84_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-38_S89_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-12_S81_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-19_S83_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-7_S80_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-35_S85_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-39_S87_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-14_S90_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-4_S78_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-37_S86_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACControl-40_S88_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-17_S82_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/AcerACDisease-6_S79_R2_clean_merged.fastq.gz \
--CPU 64 --max_memory 400G --output trinity_out_dir_AllAcerSamples_Lane1-8
```
 
### Montastrea cavernosa
Estimate total file size to estimate runtime.
```
du -ch Mcav*R1* # 53G
du -ch Mcav*R2* # 53G
```

```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/Trinity --seqType fq --left \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-24_S54_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-1_S52_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-33_S58_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-26_S46_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-18_S44_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-30_S56_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-3_S50_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-16_S48_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-32_S57_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-6_S47_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-38_S60_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-21_S45_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-3_S49_R1_clean_merged.fastq.gz  \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-36_S59_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-12_S51_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-28_S55_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-9_S43_R1_clean_merged.fastq.gz  \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-17_S53_R1_clean_merged.fastq.gz \
--right /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-24_S54_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-1_S52_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-33_S58_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-26_S46_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-18_S44_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-30_S56_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-3_S50_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-16_S48_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-32_S57_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-6_S47_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Control-38_S60_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-21_S45_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-3_S49_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Control-36_S59_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC11Disease-12_S51_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Control-28_S55_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC1Disease-9_S43_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/McavMC36Disease-17_S53_R2_clean_merged.fastq.gz \
 --CPU 64 --max_memory 400G --output trinity_out_dir_AllMcavSamples_Lane1-8 
```

### Orbicella faveolata
Estimate total file size to estimate runtime.

```
du -ch Ofav*R1* # 94G
du -ch Ofav*R2* # 93G
```
> Removing these samples from the trinity assembly because the server was running out of ram. These samples were selected for removal because they were a third replicate of a genotype. 
**R1**
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav17-7Control-40_S18_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav17-7Disease-9_S23_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofavf15Control-31_S35_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavF59Disease-19_S28_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS313Control-35_S36_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS313Disease-7_S32_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS326Control-39_S40_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS326Disease-5_S38_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Control-25_S41_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Control-34_S37_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Disease-18_S34_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Disease-3_S39_R1_clean_merged.fastq.gz \

**R2**
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav17-7Control-40_S18_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav17-7Disease-9_S23_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofavf15Control-31_S35_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavF15Disease-7_S26_R2_clean_merged.fastq.gz \
 /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavF59Disease-19_S28_R2_clean_merged.fastq.gz \
 /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS313Control-35_S36_R2_clean_merged.fastq.gz \
 /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS313Disease-7_S32_R2_clean_merged.fastq.gz \
 /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS326Control-39_S40_R2_clean_merged.fastq.gz \
 /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/OfavS326Disease-5_S38_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Control-25_S41_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Control-34_S37_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Disease-18_S34_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS326Disease-3_S39_R2_clean_merged.fastq.gz \

```
nohup /home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/Trinity --seqType fq --left \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Control-37_S16_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Control-38_S17_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Disease-19_S24_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Disease-1_S22_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Control-22_S19_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Control-26_S92_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Disease-14_S27_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Disease-2_S25_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Control-28_S20_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Control-34_S42_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Disease-12_S29_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Disease-16_S30_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Control-23_S21_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Control-29_S91_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Disease-16_S33_R1_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Disease-1_S31_R1_clean_merged.fastq.gz \
--right /home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Control-37_S16_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Control-38_S17_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Disease-19_S24_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/Ofav17-7Disease-1_S22_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Control-22_S19_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Control-26_S92_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Disease-14_S27_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF15Disease-2_S25_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Control-28_S20_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Control-34_S42_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Disease-12_S29_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavF59Disease-16_S30_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Control-23_S21_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Control-29_S91_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Disease-16_S33_R2_clean_merged.fastq.gz \
/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav/OfavS313Disease-1_S31_R2_clean_merged.fastq.gz \
--CPU 64 --max_memory 400G --output Ofav_trinity_output
 ```
 </details>
 
<details>

<summary>BBSplit</summary>
# BBSplit

1/8/24: After digging into bbsplit methods (which involved a lot of optimization) here is what was performed. 

Reliable Algal symbiont references were concatendated. 
Algal symbiont databases were made based on concatenated clade specific references.
Trinity metatranscriptomes were blastn against symbiont databases.
mapped transcripts were quality filtered.
the sequences of QC passed transcripts were extracted via cdbyank from the metatranscritpomes to make clade only references. 
These clade only references were input for bbsplit. 

# End 1/8/24 Notes to self.

> Aligns metatranscriptomes to multiple references simultaneously.

> REFERENCE Transcriptomes: (From Kelsey Beavers) The publicly available data used in this study include the transcriptomes for Symbiodinium CassKB8 (transcriptome assembly: http://medinalab.org/zoox/, accession number PRJNA80085), Breviolum minutum (transcriptome assembly: http://zoox.reefgenomics.org/download/, accession number PRJNA274852), Cladocopium goreaui (transcriptome assembly: http://ssid.reefgenomics.org/download/, accession number PRJNA307543) and Durusdinium trenchii (transcriptome assembly: https://datadryad.org/stash/dataset/doi:10.5061/dryad.12j173m, accession number PRJNA508937), as well as the genomes for M. cavernosa (genome assembly: https://matzlab.weebly.com/data-code.html, accession number PRJNA679067) and O. faveolata (genome assembly: https://www.ncbi.nlm.nih.gov/genome/13173?genome_assembly_id=311351, accession number PRJNA381078). The Master Coral database used in this study is available in a public Zenodo repository https://doi.org/10.5281/zenodo.783898080.
 
 ### Algal Symbiont
 > These references were selected because they were the refseq or highest quality reference available for the organism. 
 **Symbiodinium (Clade A):** https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001939145.1/ - GCA_001939145.1_ASM193914v1_genomic.fna.gz   
 **Breviolum minutum:** https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA274852
 **Cladocopium goreaui:** https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA307543
 **Durusdinium trenchii:** https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA508937
 
Symbiodinium file contains "@" instead of ">" at the beginning of the transcript ID. so we need to replace the @ with >
***example:***
```
sed 's/@/>/g' input.fa > output.fa
```
***modified:***
```
sed 's/@/>/g' symbiodinium_PRJNA80085.fasta  > symbiodinium_PRJNA80085.fasta 
```

### CORAL ONLY

Master Coral Database was made using the highest quality available species genomes.
 
 
**Orbicella faveolata**  https://zenodo.org/records/10151798 - Orbicella_faveolata_gen_17.scaffolds.fa # Additional genomes - https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1920453
**Montastrea cavernosa** https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz?file_subpath=%2FMcav_genome%2FMcavernosa_July2018.fasta
**Acropora cervicornis** https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/
 		Acropora cervicornis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/ - GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz 
 		Acropora cervicornis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/043/185/GCA_037043185.1_Acerv_M5/ - GCA_037043185.1_Acerv_M5_genomic.fna.gz
**Porites astreoides**
 		Porites astreoides - https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR19144705&display=download - SRR19144705.fasta.gz
 		Porites lutea - https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_958299805.1/
 		Porites australiensis - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/179/025/GCA_022179025.1_Paus_1.0/ - GCA_022179025.1_Paus_1.0_genomic.fna.gz 
 		Porites lobata - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/035/GCA_942486035.1_PLOB_v1/ - GCA_942486035.1_PLOB_v1_genomic.fna.gz  
 		Porites evermanni - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/025/GCA_942486025.1_PEVE_v1/ - GCA_942486025.1_PEVE_v1_genomic.fna.gz   
 		Porites rus - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/290/455/GCA_900290455.1_Prus/ - GCA_900290455.1_Prus_genomic.fna.gz  	
 		
To download the genome, .zip files were not downloading properly from NCBI. So, in the NCBI link for each genome, click the three vertical dots next to: "Submitted GenBank assembly
GCA_958299805.1" and select "FTP", then select the .fna.gz file to download the genome. 

Files moved from local computer to /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral on holocron server.
 
 ### Concatenate genomes to create a Master Coral Reference
 ```
 cat GCA_022179025.1_Paus_1.0_genomic.fna \
    GCA_032359415.1_NEU_Acer_K2_genomic.fna \
    GCA_037043185.1_Acerv_M5_genomic.fna \
    GCA_900290455.1_Prus_genomic.fna \
    GCA_942486025.1_PEVE_v1_genomic.fna \
    GCA_942486035.1_PLOB_v1_genomic.fna \
    GCA_958299805.1_jaPorLute2.1_alternate_haplotype_genomic.fna \
    Mcavernosa_July2018.fasta \
    Orbicella_faveolata_gen_17.scaffolds.fa \
    > concatenated_genomes.fasta
```

Renaming to MasterCoral_db
```
 /home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in concatenated_genomes.fasta -parse_seqids -dbtype nucl -out MasterCoral_db
```

Gzipping genomes to reduce space.
```
gzip GCA_022179025.1_Paus_1.0_genomic.fna
gzip GCA_032359415.1_NEU_Acer_K2_genomic.fna
gzip GCA_037043185.1_Acerv_M5_genomic.fna
gzip GCA_900290455.1_Prus_genomic.fna
gzip GCA_942486025.1_PEVE_v1_genomic.fna
gzip GCA_942486035.1_PLOB_v1_genomic.fna
gzip GCA_958299805.1_jaPorLute2.1_alternate_haplotype_genomic.fna
gzip Mcavernosa_July2018.fasta
gzip Orbicella_faveolata_gen_17.scaffolds.fa
```
 
  
Transfering references from local computer to server, I have added the NCBI accession ID to the file so the source is traceable once the reference is on the server. 
```
 scp SRR278715.fasta.gz nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/symbiodinium_PRJNA80085.fasta.gz
 scp Symbiodinium_minutum.tar.gz nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/breviolum_PRJNA274852.tar.gz
 scp CladeC_Symbiodinium_transcriptome.tar.gz nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/cladicopium_PRJNA307543.tar.gz
 scp Dtrenchii_rnaseq_assembly_v1.0.fasta nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/durusdinium_PRJNA508937.fasta
 scp MasterCoral.fasta.gz nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/MasterCoral.fasta.gz
```

###BACTERIA References:###
I Followed the instructions in step 2 [here](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#:~:text=To%20use%20the%20download%20service,button%20to%20start%20the%20download)
Essentially, this downloaded 38,177 bacteria genomes (47.5 gb) that fit this criteria Bacteria, RefSeq, Complete Genomes. 
Ideally this can be used to blast our metatranscriptomes against to identify bacteria transcripts. 
# A concise [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965038/) on the different available categories of bacteria genomes on NCBI: 
# [Overview of DIAMOND and MEGAN](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.59)

 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
> Runtime: noon April 4th
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
Concatenate all the bacteria genomes into one multi-fasta. then use makeblastdb to make that multi-fasta a reference database.
```
scp ./genome_assemblies_genome_fasta.tar nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/bacteria_reference/
```
> Runtime: 40 min. 
# decompress and concatenate with zcat
```
zcat *.fna.gz > concatenated_bacteria_genomes.fna
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in concatenated_bacteria_genomes.fna -parse_seqids -dbtype nucl -out MasterBacteria_db

/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/BacteriaOnly.Trinity.txt -num_threads 20
```


To count the number of transcripts in a fasta file
```
grep -c ">" filename.fasta
```
</details>

<details>

<summary>Longest Isoform</summary>

Create a database with blast and pull out coral only transcripts. 
First identify the longest transcript isoform in the test trinity run "trinity_out_dir_OneSample_Lanes1-8"
**Example:**
```
trinityrnaseq-Trinity-v2.5.1/util/misc/get_longest_isoform_seq_per_trinity_gene.pl raw_transctipome.fa > processed_transcriptome.fa
```
**Modified with our data:**
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_OneSample_Lanes1-8/trinity_out_dir.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_OneSample_Lanes1-8/trinity_out_dir.LongestIsoform.Trinity.fasta
```
> Runtime ~2 min. 

### Acer
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir_AllAcerSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta
```
> Runtime: 2 min
### Past
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllPastSamples_Lane1-8/trinity_out_dir_AllPastSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllPastSamples_Lane1-8/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta
```

### Ofav
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta
```

### Mcav
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllMcavSamples_Lane1-8/trinity_out_dir_AllMcavSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllMcavSamples_Lane1-8/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta
```
</details>

<details>

<summary> Make Coral Only Database</summary>
# Make Coral Only Database

Make the master coral fasta database:
**Example**
```
ncbi-blast-2.2.27+/bin/makeblastdb -in MasterCoral.fasta -parse_seqids -dbtype prot -out MasterCoral_db
```
**Modified with our data**
```
gunzip MasterCoral.fasta.gz 
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in MasterCoral.fasta -parse_seqids -dbtype prot -out MasterCoral_db
```

I am also going to make a nulceotide MasterCoral_db because the protein one takes more time for blastx. 
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in MasterCoral.fasta -parse_seqids -dbtype nucl -out MasterCoral__nucl_db
```
> Runtime ~2 min. 

**Note** my dbtype is prot for protein. This means it is producing a protein master coral database. so when I do my blast below, I have a nucleotide query (trinity_out_dir.LongestIsoform.Trinity.fasta) mapping to a protein database so I need to use blastx. If I made the master coral db a dbtype "nucl" for nucleotide than I could use a blastn blast. 

Here is a good [visual](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/) on blast types.

```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastx -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_OneSample_Lanes1-8/trinity_out_dir.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_OneSample_Lanes1-8/trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt
```
</details>

<details>

<summary>Blastn</summary>

### Acer
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastx -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt
```
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt -num_threads 20
```

### Past 
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20
```
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/LongestIsoform.BacteriaOnly.Trinity.txt -num_threads 20
```

### Ofav
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.CoralOnly.Trinity.txt -num_threads 20
```
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.BacteriaOnly.Trinity.txt -num_threads 20
```

### Mcav
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20
```
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/BacteriaOnly.Trinity.txt -num_threads 20
```
</details>

<details>
<summary>Filter</summary>
Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
###Acer - Coral Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt > Acer_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Acer_contigs_percent_95.txt > Acer_contigs_percent_95_bp_150.txt
```
###Acer - Bacteria Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' BacteriaOnly.Trinity.txt > Acer_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Acer_Bacteria_contigs_percent_95.txt > Acer_Bacteria_contigs_percent_95_bp_150.txt
```

###Past - Coral Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt > Past_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Past_contigs_percent_95.txt > Past_contigs_percent_95_bp_150.txt
```

###Past - Bacteria Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' LongestIsoform.BacteriaOnly.Trinity.txt > Past_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Past_Bacteria_contigs_percent_95.txt > Past_Bacteria_contigs_percent_95_bp_150.txt
```

###Mcav - Coral Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt > Mcav_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Mcav_contigs_percent_95.txt > Mcav_contigs_percent_95_bp_150.txt
```

###Mcav - Bacteria Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' LongestIsoform.BacteriaOnly.Trinity.txt > Mcav_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Mcav_Bacteria_contigs_percent_95.txt > Mcav_Bacteria_contigs_percent_95_bp_150.txt
```

###Ofav - Coral Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' Ofav_trinity_output.LongestIsoform.CoralOnly.Trinity.txt > Ofav_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Ofav_contigs_percent_95.txt > Ofav_contigs_percent_95_bp_150.txt
```

###Ofav - Bacteria Only
```
awk '{if ($3 > 95) print $1,$2,$4 }' Ofav_trinity_output.LongestIsoform.BacteriaOnly.Trinity.txt > Ofav_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Ofav_Bacteria_contigs_percent_95.txt > Ofav_Bacteria_contigs_percent_95_bp_150.txt
```
</details>

<details>
<summary>Transdecoder</summary>

> Transdecoder first appeared as a straight-forward "convert nucleotides to protein prediction" intermediate step. I quickly realized a significant overlooked issue with this that I have addressed in my methods below.

> The Problem: Orthofinder relies on the predicted proteins from Transdecoder. To predict proteins, transdecoder identifies Open Reading Frames "ORFs". These are the nucleotides that create the start codon, so a logical way to predict the start of a protein. The problem is that there can be multiple, possibly hundreds of ORFs per gene. So, Orthofinder then using these multiple ORFs per gene to identify orthogroups, single copy orthologs, etc. These multiple ORF predicted orthogroups become a problem when its time to merge that information with count data produced from Salmon. Salmon uses the reference transcriptome that we used from BBsplit to count the frequency of those genes. So these are gene-level counts, NOT multiple ORF within a gene-level counts. So what do you do when you have lets say 750 ORFs, with unique functions, molecular evolution and need to assign a count at the gene-level? The answer is it doesnt work cleanly and creates a technical issue. There are various ways people address this and I played with most of them. You can pick the first ORF, the longest ORF, the ORF with the most homology. When I picked the first ORF, the problem is you have to just delete the other ORFs and that feels so messy and arbitrary. When I picked the longest ORF, I observed very strong species bias. This is probably because across mutliple species, some may evolutionarily have duplication events in their genes that will make them longer, or their references are higher quality leading to better alignment because its a better studied speices or their genome happened to be produced more recently with better equipment, or something about that species does better during the benchwork stage and extracts better than the others. The first function of transdecoder: transdecoder.LonOrfs, states it selects the longest ORF, but when you read the details of this function, it does not actually ONLY retain the longest ORFs and has exceptions that will have it keep multiple ORFs. I wrote a python script to force it to retain the longest ORF by identifying the transcript/gene names that were the same, and then keeping the one that was the longest length. These details can be found in the transcript titles. Once I did this it gave me the species bias I described and was not purusued. I do not include the script to force the retention of longest ORF. All of these are technical artificats and not a bioliogical factor that justifies the longest ORF species bias approach. The final approach is a ORF homology-based approach. During the Transdecoder.Predict function, you can specify "single_best_only". This function will select only one representative ORF per transcript based on homology. During the prediction process, TransDecoder evaluates ORFs based on coding potential, length, and homology to known protein domains or sequences.
ORFs are ranked by these scores, with higher scores indicating a stronger likelihood of representing true protein-coding sequences as the "single-best" ORF. I found this homolgy based approach to represent the biological nature of the data while the others introduced a technical bias. 

> When not to select one ORF: Selecting one ORF is ideal for those who are interested in gene or transcript-level analysis and do not care about isoform or ORF level analysis. Much of multi-species comparisons is a tradeoff between complexity and comparability. If you want to do gene-level comparisons, there is likely to be more conservation of those genes and they are less species-specific and comparable. If you want to do ORF-level comparisons, you will retain the complexity of each species, but that complexity may be species-specific and less evolutionarity conserved. Because we are comparing four species in this study that expand potentially 100+ Million years of evolutionary divergence, I am leaning into comparability. But if I was compariing Orbicella faveolata to Orbicella franksi or even Orbicella faveolata to Montastraea cavernosa (15 million year divergence), I may lean into retaining that complexity to explore the molecular evolution of these relatively closely related species. 

### *Retaining Non-Protein-Coding Regions*

> BONUS. From this deep dive into protein-prediction, I realized the default transdecoder only predicts protein-coding regions and not non-protein-coding regions. Why should we care about non-protein-coding regions? Non-protein-coding regions are highly valuable in multi-species RNA expression comparison studies, particularly in the context of disease response, because they often regulate gene expression, mediate cellular signaling, and adaptively respond to environmental or pathological stimuli. Because our questions evolve around what is SCTLD and why are some species resistnat, these non-protein-coding regions are worth retaining and this pipeline integrates how to do this.
To briefly explain, in adition to protein-coding protein prediction with transdecoder, we also perform pfam and a blastp search. This gives us three "datasets" of semi-overlapping results that we feed into TransDecoder.Predic and have it retain the "single_best_only" per transcript. When you just retain protein-coding regions, you are inadvertenly discarding about ~60% of the transcriptome. So implementing pfam and blastp allows you to maximize your dataset. 

### Verify single ORF per transcript with a python script:
> This command is helpful at verifying your output is only one ORF per transcript.

```
nano
```
copy and paiste this script into nano:
```
#!/usr/bin/env python3

import sys
import re
from collections import defaultdict

# Check if a FASTA file argument is provided
if len(sys.argv) != 2:
    print("Usage: ./check_single_orf.py <fasta_file>")
    sys.exit(1)

# Get the FASTA file from the command-line arguments
fasta_file = sys.argv[1]

# Dictionary to store ORF counts per transcript
transcript_orf_counts = defaultdict(list)

# Regular expression to extract the transcript ID before the ORF identifier (.pX)
orf_pattern = re.compile(r"^(>)(\S+?)\.p\d+\b")  # Matches '>TRINITY_ID.pX'

# Parse the FASTA file
try:
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                match = orf_pattern.search(line)
                if match:
                    transcript_id = match.group(2)  # Extract the base transcript ID
                    transcript_orf_counts[transcript_id].append(line.strip())
except FileNotFoundError:
    print(f"Error: File '{fasta_file}' not found.")
    sys.exit(1)

# Check for transcripts with multiple ORFs
multiple_orfs = {tid: orfs for tid, orfs in transcript_orf_counts.items() if len(orfs) > 1}

# Report results
if multiple_orfs:
    print("Transcripts with multiple ORFs detected:")
    for transcript, orfs in multiple_orfs.items():
        print(f"\nTranscript ID: {transcript}")
        for orf in orfs:
            print(f"  {orf}")
else:
    print("Verification passed: Each transcript has only a single ORF.")
```


```
chmod +x check_single_orf.py
```

# Bacteria :microbe:
## Acer - Bacteria 
```
mkdir Bacteria_transdecoder_AllORFs
cd Bacteria_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Acer_Bacteria_only_transcriptome.fa
```
Pfam using Hmmer
```
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

```
Blastp
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Bacteria_transdecoder_AllORFsAcer_Bacteria_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6
```
Transdecoder.Predict with Transdecoder.LongORFs, Pfam, and blastp results and retaining --single_best_only ORF per transcript.
```
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../../Acer_Bacteria_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only
```
Combine similar transcripts into one using CD-HIT
```
mv Acer_Bacteria_only_transcriptome.fa.transdecoder.pep Acer_Bacteria_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Acer_Bacteria_only_transcriptome_transdecoder.fa -o Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
```

*Verify single ORF per transcript:*
> Update "fasta.file" to be your output file from transdecoder. So for the Acer-Bacteria example below, I would replace "file.fasta" with "Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa"
```
./check_single_orf.py file.fasta
```


## Mcav - Bacteria
```
mkdir Bacteria_transdecoder_AllORFs
cd Bacteria_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Mcav_Bacteria_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Bacteria_transdecoder_AllORFs/Mcav_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Bacteria_transdecoder_AllORFs/Mcav_Bacteria_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Mcav_Bacteria_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Mcav_Bacteria_only_transcriptome.fa.transdecoder.pep Mcav_Bacteria_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Mcav_Bacteria_only_transcriptome_transdecoder.fa -o Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
```
## Ofav - Bacteria
```
mkdir Bacteria_transdecoder_AllORFs
cd Bacteria_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Ofav_Bacteria_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Bacteria_transdecoder_AllORFs/Ofav_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Bacteria_transdecoder_AllORFs/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Ofav_Bacteria_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Ofav_Bacteria_only_transcriptome.fa.transdecoder.pep Ofav_Bacteria_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Ofav_Bacteria_only_transcriptome_transdecoder.fa -o Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
```
## Past - Bacteria
```
mkdir Bacteria_transdecoder_AllORFs
cd Bacteria_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Past_Bacteria_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Bacteria_transdecoder_AllORFs/Past_Bacteria_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Bacteria_transdecoder_AllORFs/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Past_Bacteria_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Past_Bacteria_only_transcriptome.fa.transdecoder.pep Past_Bacteria_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Past_Bacteria_only_transcriptome_transdecoder.fa -o Past_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Past_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
```

# Host :house:

## Acer - Host
```
/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer
mkdir Host_transdecoder_AllORFs
cd Host_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Acer_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Host_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Host_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Host_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

mv blastp.outfmt6 ../


# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Acer_coral_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Acer_coral_only_transcriptome.fa.transdecoder.pep Acer_coral_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Acer_coral_only_transcriptome_transdecoder.fa -o Acer_Host_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Acer_Host_reference_proteome_AllORF_SingleBestOnly.fa
```

## Mcav - Host
```
mkdir Host_transdecoder_AllORFs
cd Host_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Mcav_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Host_transdecoder_AllORFs/Mcav_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Host_transdecoder_AllORFs/Mcav_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Host_transdecoder_AllORFs/Mcav_coral_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

mv blastp.outfmt6 ../
# Rename
mv Mcav_coral_only_transcriptome.fa.transdecoder_dir/ Mcav_coral_only_transcriptome.fa.transdecoder_dir

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Mcav_coral_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Mcav_coral_only_transcriptome.fa.transdecoder.pep Mcav_coral_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Mcav_coral_only_transcriptome_transdecoder.fa -o Mcav_coral_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Mcav_coral_reference_proteome_AllORF_SingleBestOnly.fa
```

## Ofav - Host
```
mkdir Host_transdecoder_AllORFs
cd Host_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Ofav_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Host_transdecoder_AllORFs/Ofav_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Host_transdecoder_AllORFs/Ofav_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Host_transdecoder_AllORFs/Ofav_coral_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Ofav_coral_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Ofav_coral_only_transcriptome.fa.transdecoder.pep Ofav_coral_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Ofav_coral_only_transcriptome_transdecoder.fa -o Ofav_Host_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Ofav_Host_reference_proteome_AllORF_SingleBestOnly.fa
```

## Past - Host
```
mkdir Host_transdecoder_AllORFs
cd Host_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Past_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Host_transdecoder_AllORFs/Past_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Host_transdecoder_AllORFs/Past_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Host_transdecoder_AllORFs/Past_coral_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Past_coral_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Past_coral_only_transcriptome.fa.transdecoder.pep Past_coral_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Past_coral_only_transcriptome_transdecoder.fa -o Past_Host_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Past_Host_reference_proteome_AllORF_SingleBestOnly.fa
```
# Algae :seedling:

## Acer - Algae
```
/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Acer_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Acer_coral_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Acer_Algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Acer_Algae_only_transcriptome.fa.transdecoder.pep Acer_Algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Acer_coral_only_transcriptome_transdecoder.fa -o Acer_Algae_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Acer_Algae_reference_proteome_AllORF_SingleBestOnly.fa
```

## Mcav - Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Mcav_Algae_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Mcav_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Mcav_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Mcav_Algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Mcav_Algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Mcav_Algae_only_transcriptome.fa.transdecoder.pep Mcav_Algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Mcav_Algae_only_transcriptome_transdecoder.fa -o Mcav_Algae_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Mcav_Algae_reference_proteome_AllORF_SingleBestOnly.fa
```

## Ofav - Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Ofav_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Algae_transdecoder_AllORFs/Ofav_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs/Ofav_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs/Ofav_Algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Ofav_Algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Ofav_Algae_only_transcriptome.fa.transdecoder.pep Ofav_Algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Ofav_Algae_only_transcriptome_transdecoder.fa -o Ofav_Algae_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Ofav_Algae_reference_proteome_AllORF_SingleBestOnly.fa
```

## Past - Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Past_coral_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Past_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Past_Algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Past_Algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Past_Algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Past_Algae_only_transcriptome.fa.transdecoder.pep Past_Algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Past_Algae_only_transcriptome_transdecoder.fa -o Past_Algae_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Past_Algae_reference_proteome_AllORF_SingleBestOnly.fa
```

# OrthoFinder

### Move all the reference_proteome.fa into a new folder specific for each holobiont compartment (Host, Algae, Bacteria, so three folders total)
Example of moving all reference_proteome.fa files into their respective newly created folder:
```
# Host
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF
mv *_Host_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/

# Algae
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF
mv *_Algae_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/

# Bacteria
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF
mv *_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/
```

Then within each folder as the active directory perform the Orthofinder Command:
> This command performs OrthoFinder, comparing the predicted proteins among the reference_proteome.fa files
```
# Host
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50

# Algae
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50

# Bacteria
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50
```

</details>




<details>
<summary>BBMAP</summary>
# BBMAP
> Reference genomes need to be indexed before bbsplit can be ran. 

```
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa build=1
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa build=2
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/symbiodinium_GCA_001939145.1.fna build=
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/breviolum_PRJNA274852.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/durusdinium_PRJNA508937.fasta
 ```
</details>

<details>
<summary>BBSplit</summary>
# BBSplit

### Applied by For Loop. This doesnt work perfectly, it runs but doesnt apply to all samples and is inconsistent so I manually ran bbsplit which is referenced below. 
 ```
#!/bin/bash

# Add directories to the PATH
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/

# Define the reference directory
DIR=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references
HDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer
SDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer

# Process each FASTQ file
find ${SDIR} -name "*_R1_clean_merged.fastq.gz" | while read FILE; do
    echo ${FILE}
    SAMP=$(basename ${FILE} _R1_clean_merged.fastq.gz)
    echo ${SAMP}
nohup /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh in1="${SDIR}/${SAMP}_R1_clean_merged.fastq.gz" in2="${SDIR}/${SAMP}_R2_clean_merged.fastq.gz" ref="${HDIR}/Acer_coral_only_transcriptome.fa,${HDIR}/Acer_Bacteria_only_transcriptome.fa,${DIR}/symbiodinium_PR$
basename="${SAMP}_%.fq.gz" refstats="${SAMP}_stats.txt" ambig=all ambig2=all
outu1="${SAMP}_bboutu_1.fq.gz" outu2="${SAMP}_bboutu_2.fq.gz" &
done
```
### Applied Manually:
> Two samples as example, this will need to be scaled up for each sample, or get the for loop adapted to your project. 

```
nohup /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh in1="/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer/AcerACControl-37_S86_R1_clean_merged.fastq.gz" in2="/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer/AcerACControl-37_S86_R2_clean_merged.fastq.gz" ref="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa,/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/symbiodinium_GCA_001939145.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/breviolum_PRJNA274852.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/durusdinium_PRJNA508937.fasta" basename="AcerACControl-37_S86_%.fq.gz" refstats="AcerACControl-37_S86_stats.txt" ambig=all ambig2=all outu1="AcerACControl-37_S86_bboutu_1.fq.gz" outu2="AcerACControl-37_S86_bboutu_2.fq.gz"
mv AcerAC* Acer_output
rm nohup.out
rm -r ref/
nohup /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh in1="/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer/AcerACDisease-14_S90_R1_clean_merged.fastq.gz" in2="/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer/AcerACDisease-14_S90_R2_clean_merged.fastq.gz" ref="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa,/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/symbiodinium_GCA_001939145.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/breviolum_PRJNA274852.fa,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta,/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/durusdinium_PRJNA508937.fasta" basename="AcerACDisease-14_S90_%.fq.gz" refstats="AcerACDisease-14_S90_stats.txt" ambig=all ambig2=all outu1="AcerACDisease-14_S90_bboutu_1.fq.gz" outu2="AcerACDisease-14_S90_bboutu_2.fq.gz"
```

</details>






<details>
<summary> BBSplit - Seperate F and R Reads</summary>
	
> Split bbsplit output into Forward and Reverse reads using BBmap

**Example:**
```
#!/bin/bash
PATH=$PATH:/opt/storage/opt_programs/bbmap/
PATH=$PATH:/opt/storage/anaconda3/bin/java
for FILE in *_host.fq.gz; do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
nohup /opt/storage/anaconda3/bin/java -ea -Xmx10g  -cp  /opt/storage/opt_programs/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=${SAMP}_1.fq out2=${SAMP}_2.fq
done
```

**Modified - host**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_Acer_coral_only_transcriptome.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
done
```
**Modified - bacteria**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_Acer_Bacteria_only_transcriptome.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
done
```
**Modified - symbiodinium_GCA_001939145**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_symbiodinium_GCA_001939145.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
done
```

**Modified - breviolum_PRJNA274852**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_breviolum_PRJNA274852.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
done
```

**Modified - davies_cladeC_feb**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_davies_cladeC_feb.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
done
```

**Modified - durusdinium_PRJNA508937**
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in ./$(ls *_durusdinium_PRJNA508937.fq.gz); do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=../Acer_output_FR/${SAMP}_1.fq out2=../Acer_output_FR/${SAMP}_2.fq
Done
```
</details>

<details>
<summary>Salmon Indexing</summary>
# Salmon
> Index and then Read Quantification

**Salmon Indexing** 
First have to build a salmon index for your transcriptome. Assume that transcripts.fa contains the set of transcripts you wish to quantify. 

### Acer
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Acer/Acer_index
```

### Mcav
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Mcav/Mcav_index
```

### Past
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Past/Past_index
```

### Ofav
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Ofav/Ofav_index
```

## Index Building for Symbiont Transcriptomes 

### Symbiodinium index - with reference genome
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/symbiodinium_GCA_001939145.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_A_index -k 23
```

### Symbiodinium index - with bbsplit output transcriptome
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/durusdinium_PRJNA508937.fasta -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_A_index -k 23
```

### Breviolium index
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/breviolum_PRJNA274852.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_B_index -k 23
```

### Cladicopium index
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_C_index -k 23
```

### Durusdinium index
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/durusdinium_PRJNA508937.fasta -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_D_index -k 23
```


### Index Building for Bacteria Transcriptomes
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/bacteria_reference/concatenated_bacteria_genomes.fna -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Bacteria_index -k 23
```

#### Salmon Example:
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Clade_A_index -l A -1 /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Acer/Acer_output_FR/AcerACDisease-7_S80_symbiodinium_GCA_001939145_1.fq -2 /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Acer/Acer_output_FR/AcerACDisease-7_S80_symbiodinium_GCA_001939145_2.fq -p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/AcerACDisease-7_S80_symbiodinium_quant
```

count number of transcripts with a count equal or greater to 1
```
awk 'NR>1 && $5 >= 1 { count++ } END { print "Number of transcripts with NumReads >= 1: ", count }' quant.sf
```
</details>
<details>
<summary>Salmon - Quantification</summary>
# Salmon - Quantification

 
### Loop for Acer
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Acer_coral_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Acer_coral_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Acer/Acer_index -l A \
				-1 ${SAMP}_Acer_coral_only_transcriptome_1.fq \
				-2 ${SAMP}_Acer_coral_only_transcriptome_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_host_quant
done
```

### Loop for Symbiodinium
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_symbiodinium_GCA_001939145_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _symbiodinium_GCA_001939145_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Clade_A_index -l A \
				-1 ${SAMP}_symbiodinium_GCA_001939145_1.fq \
				-2 ${SAMP}_symbiodinium_GCA_001939145_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_Clade_A_quant
done
```

### Loop for Breviolium
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_breviolum_PRJNA274852_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _breviolum_PRJNA274852_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Clade_B_index -l A \
				-1 ${SAMP}_breviolum_PRJNA274852_1.fq \
				-2 ${SAMP}_breviolum_PRJNA274852_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_Clade_B_quant
done
```

### Loop for Cladicopium
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_davies_cladeC_feb_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _davies_cladeC_feb_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Clade_B_index -l A \
				-1 ${SAMP}_davies_cladeC_feb_1.fq \
				-2 ${SAMP}_davies_cladeC_feb_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_Clade_C_quant
done
```

### Loop for Durusdinium
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_durusdinium_PRJNA508937_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _durusdinium_PRJNA508937_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Clade_B_index -l A \
				-1 ${SAMP}_durusdinium_PRJNA508937_1.fq \
				-2 ${SAMP}_durusdinium_PRJNA508937_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_Clade_D_quant
done
```

### Loop for Bacteria
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Acer_Bacteria_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Acer_Bacteria_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Clade_B_index -l A \
				-1 ${SAMP}_Acer_Bacteria_only_transcriptome_1.fq \
				-2 ${SAMP}_Acer_Bacteria_only_transcriptome_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/${SAMP}_Bacteria_quant
done
```
### Repeat for other Coral

</details>

<details>
<summary>Orthofinder</summary>
# Orthofinder: Obtaining single-copy orthologs

### Need to make reference proteomes for the Algae.

## Clade A
```
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t symbiodinium_GCA_001939145.fa
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t symbiodinium_GCA_001939145.fa
mv symbiodinium_GCA_001939145.fa.transdecoder.pep symbiodinium_GCA_001939145_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i symbiodinium_GCA_001939145_transdecoder.fa -o Clade_A_reference_proteome.fa
```

## Clade B
```
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t breviolum_PRJNA274852.fa
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t breviolum_PRJNA274852.fa
Mv breviolum_PRJNA274852.fa.transdecoder.pep Clade_B_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_B_transdecoder.fa -o Clade_B_reference_proteome.fa
```

## Clade C
```
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t CladeC_Symbiodinium_transcriptome/davies_cladeC_feb.fasta
mv davies_cladeC_feb.fasta.transdecoder.pep Clade_C_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_C_transdecoder.fa -o Clade_C_reference_proteome.fa
```

## Clade D
```
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t durusdinium_PRJNA508937.fasta
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t durusdinium_PRJNA508937.fasta
mv durusdinium_PRJNA508937.fasta.transdecoder.pep Clade_D_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_D_transdecoder.fa -o Clade_D_reference_proteome.fa
```

### Move all proteome files to proteome folder. Then run Orthofinder as follows:

```
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 20
```
> This will generate a directory called "OrthoFinder/Results\_[date]". Move the "Orthogroups" and the "Comparative_Genomics_Statistics" subdirectories onto your local machine.

</details>

<details>

<summary> Annotating the Orthologs</summary>
# Annotating the orthologs
From the /Orthogroups directory, use "Orthogroups_SingleCopyOrthologues.txt" and "Orthogroups.tsv" to grab the sequence names of the single-copy orthologs from each species.
```{r}
setwd("~/OneDrive - University of Texas at Arlington/Dissertation/SCTLD/Files_for_R/January_2022/SCTLD/host/all/Orthofinder_cdhit/Orthogroups/")
Orthogroups_SingleCopyOrthologues <- read.csv("Orthogroups_SingleCopyOrthologues.txt",header = FALSE)
colnames(Orthogroups_SingleCopyOrthologues) <- c("Orthogroup")
Orthogroups <- read.table(file = 'Orthogroups.tsv', sep = '\t', header = TRUE)
Orthogroup_Sequences <- merge(Orthogroups_SingleCopyOrthologues,Orthogroups,by="Orthogroup")
write.csv(Orthogroup_Sequences, file="Orthogroup_Sequences.csv",quote = FALSE,row.names = FALSE)

# We also grab the sequence names from P. astreoides (the species with highest quality assembly) and will use these for annotation 
names(Orthogroup_Sequences)
mcav_orthologs <- Orthogroup_Sequences[c(1,3)]
names(mcav_orthologs)
write.csv(mcav_orthologs, file="mcav_SC_ortholog_sequences.csv",quote = FALSE,row.names = FALSE)
```

Make a new directory in "OrthoFinder/" called "Annotating_Orthogroups" and move "mcav_SC_ortholog_sequences.csv" here along with the M.cavernosa proteome from the previous step. We will use [cdbfasta](https://github.com/gpertea/cdbfasta) to pull the single-copy ortholog sequences out of the M.cavernosa reference proteome.

```{linux, eval=FALSE}
DIR=/opt/storage/storage/SCTLD/host_orthofinder/OrthoFinder/OrthoFinder/Annotating_Orthogroups

# Make an index of the P.astreoides reference proteome
cdbfasta/cdbfasta mcav_ref_proteome.fa

cat mcav_SC_ortholog_sequence_names.txt | cdbfasta/cdbyank mcav_ref_proteome.fa.cidx > mcav_sc_orthologs.fa

# Annotate with BLASTp
ncbi-blast-2.2.27+/bin/blastp -query mcav_sc_orthologs.fa -db uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out mcav_sc_orthologs_annotated.txt
```



Save "mcav_sc_orthologs_annotated.txt" and "mcav_sc_orthologs.fa" to your computer. We will use these files in our Coral EVE analysis. 
This process is analagous in the Symbiodinaiceae genera using the Durusdinium reference proteome as the reference for annotation.

</details>
