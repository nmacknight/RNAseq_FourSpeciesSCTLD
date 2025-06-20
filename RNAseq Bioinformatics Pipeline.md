# RNAseq Bioinformatics Pipeline 
🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling: :microbe: 🪸 :seedling:


**Project:** ***Acropora cervicornis***, ***Porites astreoides***, ***Montastraea cavernosa***, and ***Orbicella faveolata*** were exposed to Stony Coral Tissue Loss Disease. Tissue samples were collected to investigate the gene expression of each component of the coral holobiont with particular consideration of the Bacteria Gene expression as this focus has never been investigated and may lead to results that identify a bacteria disease causing agent and mechanistic etiology. 

**Goal:** Process sequence reads to build Metatranscriptomes, isolate holobiont components (Host 🪸, Algal Symbiont :seedling:, Bacteria :microbe:), identify homologous transcripts, and annotate those orthogroups for statistical analysis in R.


> [!TIP]
> Click the arrows to expand each section




<details>

<summary>Downloading Software</summary>

# Downloading Software

> Notes for you, another user: I had user permission limitations due to this being conducted on a government server and did not and could not be granted admin privileges or the ability to use sudo. So, while this is how I installed and compiled the necessary software, if you have admin priviliges on your server, you are likely to run into less errors, hopefully. It also means the exact code to install the software may differ for you.

> Notes for myself: I am installing this software on a server which required me to install basic software into my user folder because it was a blank slate. I did find a parent folder that had a lot of software, but there were inconsistencies on whether I could execute using those pre-installed software tools and certainly could not write new software into that parent folder. Also, some software, I believe Trinity and Samtools were successfully compiled either by ethernet connection in the aoml building like how I mentioned below, and/or by being in a conda environment and then performing the compiling or running the code.


**Anaconda**
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

**FASTP**
> Fastp is a tool designed to provide fast all-in-one preprocessing for FastQ files.
```
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```

**TRINITY**
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

**BLAST**
> Extract Coral Only (or Symbiont Only or Bacteria Only) Transcripts
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
```

**cd-hit**
> Create BLAST Database
```
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
tar -xvzf cd-hit-v4.8.1-2019-0228.tar.gz
make
cd cd-hit-auxtools
make
```

**cdbfasta**
> CDBFASTA is a tool used for indexing and querying large sequence sets in the FASTA format efficiently. It creates a hash index of sequences, allowing for fast retrieval and comparison of sequences based on their headers or identifiers.
```
git clone https://github.com/gpertea/cdbfasta.git
cd cdbfasta
make
```
**UNIPROT**
> This is the database that will be used to annotate sequences and assign an Entry ID and an E-value number which is a score on how well the sequence matches the database sequence. Low E values represent high sequence matching, so you want low e-values. As a standard bench mark, E-values with less than e^-05 or 0.00001 are retained and anything greater than that is filtered out.
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```

**BOOST**
> This may have not been utilized. 
```
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
```

**TOPHAT**
```
git clone https://github.com/infphilo/tophat.git
./bootstrap
./configure
```

**DB_File**
```
cpan
install DB_File
exit
```

**SAMTOOLS**
```
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar -xjf samtools-1.14.tar.bz2
./configure --without-curses
make
```
> ERROR : Permission denied. Probably fixed by trying once connected to aoml ethernet


**BBMAP**
```
wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz
tar -xvzf BBMap_38.90.tar.gz
```

**SALMON**
```
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/Salmon-1.5.2_linux_x86_64.tar.gz
tar -xvzf Salmon-1.5.2_linux_x86_64.tar.gz
```

**OrthoFinder**
```
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder_source.tar.gz
tar xzf OrthoFinder_source.tar.gz
python OrthoFinder_source/orthofinder.py -h # to verify installation
```

**Transdecoder**
```
wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz
tar -xvzf TransDecoder-v5.7.1.tar.gz
```

**Singularity**
```
wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/archived/trinityrnaseq.v2.15.0.simg
tar -xzf trinityrnaseq.v2.15.0.simg
cd singularity
./mconfig
make -C builddir
```

**Libevent**
```
wget http://monkey.org/~provos/libevent-1.4.14b-stable.tar.gz
tar xzf libevent-1.4.14b-stable.tar.gz
cd libevent-1.4.14b-stable
./configure --prefix=/opt/libevent
make
make install #permission denied
```

**ncurses**
```
wget http://invisible-mirror.net/archives/ncurses/current/ncurses.tar.gz
tar xzf ncurses.tar.gz
cd ncurses-6.4-20240224
./configure
make
```

**tmux**
```
wget https://github.com/tmux/tmux/releases/download/3.4/tmux-3.4.tar.gz
tar xvf tmux-3.4.tar.gz
cd tmux-3.4/
./configure # needs libevent. Failed. 
make
```

**BLAST**
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.15.0+-x64-linux.tar.gz
cd ncbi-blast-2.15.0+
export PATH=$PATH:/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin
```
**BUSCO**
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

You can occassionally check on how many files are completed (after opening another terminal window) with:
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

<summary>Trinity - Metatranscriptome Assembly</summary>
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

Acer
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir_AllAcerSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta
```
> Runtime: 2 min

Past
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllPastSamples_Lane1-8/trinity_out_dir_AllPastSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllPastSamples_Lane1-8/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta
```

Ofav
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta
```

Mcav
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllMcavSamples_Lane1-8/trinity_out_dir_AllMcavSamples_Lane1-8.Trinity.fasta > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllMcavSamples_Lane1-8/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta
```
</details>

<details>

<summary> Make Coral Host Database</summary>
# Make Coral Only Database

> The convenient part about the host reference is that visually we know what the coral was, so there isnt a necessity to have a species-specific reference be our only reference and instead we can combine coral species references to capture a greater diversity of coral host transcripts within our metatranscriptomes.

.zip files were not downloading properly from NCBI. So, in the NCBI link for each genome, click the three vertical dots next to: Submitted GenBank assembly
GCA_958299805.1 and select FTP, then select the .fna.gz file to just download the genome. 

Files moved to /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral on holocron server.

Below are the locations and names of the curated coral host references. 
 
Orbicella faveolata
```
# Orbicella_faveolata_gen_17.scaffolds.fa 
 https://zenodo.org/records/10151798

 https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1920453
```
 
Montastrea cavernosa 
```
https://www.dropbox.com/s/yfqefzntt896xfz/Mcavernosa_genome.tgz?file_subpath=%2FMcav_genome%2FMcavernosa_July2018.fasta
```

Acropora cervicornis 
```
https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/

# GCA_032359415.1_NEU_Acer_K2_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/

# GCA_037043185.1_Acerv_M5_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/043/185/GCA_037043185.1_Acerv_M5/
```

Porites astreoides 
```
Porites astreoides - https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR19144705&display=download - SRR19144705.fasta.gz
Porites lutea - https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_958299805.1/
Porites australiensis - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/022/179/025/GCA_022179025.1_Paus_1.0/ - GCA_022179025.1_Paus_1.0_genomic.fna.gz 
Porites lobata - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/035/GCA_942486035.1_PLOB_v1/ - GCA_942486035.1_PLOB_v1_genomic.fna.gz  
Porites evermanni - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/025/GCA_942486025.1_PEVE_v1/ - GCA_942486025.1_PEVE_v1_genomic.fna.gz   
Porites rus - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/290/455/GCA_900290455.1_Prus/ - GCA_900290455.1_Prus_genomic.fna.gz  	
```

 
Concatenate genomes
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


```
 /home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in concatenated_genomes.fasta -parse_seqids -dbtype nucl -out MasterCoral_db
```
Gzipping them again to reduce space.
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

> Runtime ~2 min. 

> **Note** my dbtype is prot for protein. This means it is producing a protein master coral database. so when I do my blast below, I have a nucleotide query (trinity_out_dir.LongestIsoform.Trinity.fasta) mapping to a protein database so I need to use blastx. If I made the master coral db a dbtype "nucl" for nucleotide than I could use a blastn blast. 

Here is a good [visual](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/) on blast types.

Acer
```
#blastn
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt -num_threads 20

#Filter
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt > Acer_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Acer_contigs_percent_95.txt > Acer_contigs_percent_95_bp_150.txt

#Index Metatranscriptome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta

#cdbyank - Extract filtered sequences
cat Acer_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Acer_coral_only_transcriptome.fa
```

Mcav
```
#blastn
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20

#Filter
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt > Mcav_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Mcav_contigs_percent_95.txt > Mcav_contigs_percent_95_bp_150.txt

#Index Metatranscriptome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta

#cdbyank - Extract filtered sequences
cat Mcav_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Mcav_coral_only_transcriptome.fa
```

Ofav
```
#Blastn
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.CoralOnly.Trinity.txt -num_threads 20

#Filter
awk '{if ($3 > 95) print $1,$2,$4 }' Ofav_trinity_output.LongestIsoform.CoralOnly.Trinity.txt > Ofav_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Ofav_contigs_percent_95.txt > Ofav_contigs_percent_95_bp_150.txt

#Index
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Ofav_trinity_output.LongestIsoform.Trinity.fasta

#cdbyank
cat Ofav_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Ofav_trinity_output.LongestIsoform.Trinity.fasta.cidx > Ofav_coral_only_transcriptome.fa
```

Past
```
#Blastn
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20

#Filter
awk '{if ($3 > 95) print $1,$2,$4 }' trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt > Past_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Past_contigs_percent_95.txt > Past_contigs_percent_95_bp_150.txt

#Index Metatranscriptome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta

#cdbyank
cat Past_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Past_coral_only_transcriptome.fa
```

</details>

<details>

<summary> Make Algal Symbiont Database</summary>

> Decompress references for accurate concatenation.

Concatenate genomes
```
 # Clade A - Symbiodinium
 cat Symbiodinium_Aranda2012.fa \
 	Symbiodinium_Shoguchi2018.fasta \
    > concatenated_Clade_A_genomes.fasta
    
# Clade B - Breviolum
 cat Breviolum_Shoguchi2013.fa \
 	Breviolum_AvilaMagana2021.fna \
    > concatenated_Clade_B_genomes.fasta

# Clade C - Cladicopium
 cat Cladocopium_sp_C92/Cladocopium_sp_C92.genome.fa \
 	Cladocopium_goreaui/Cladocopium_goreaui.genome.fa \
    > concatenated_Clade_C_genomes.fasta
    
# Clade D - Durusdinium
 cat Durusdinium_Shoguchi2021.fa \
 	durusdinium_PRJNA508937.fasta \
    > concatenated_Clade_D_genomes.fasta
```
Make Clade Database
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb \
    -in concatenated_Clade_A_genomes.fasta \
    -parse_seqids -dbtype nucl -out Clade_A_db

/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb \
    -in concatenated_Clade_B_genomes.fasta \
    -parse_seqids -dbtype nucl -out Clade_B_db
    
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb \
    -in concatenated_Clade_C_genomes.fasta \
    -parse_seqids -dbtype nucl -out Clade_C_db
    
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb \
    -in concatenated_Clade_D_genomes_unique.fasta \
    -parse_seqids -dbtype nucl -out Clade_D_db

```
Confirm accurate concatentation by counting the number of transcripts in the individual files add up to the new concatenated reference with:
```
grep -c "^>" Symbiodinium_Shoguchi2018.fasta
# 164631
grep -c "^>" Symbiodinium_Aranda2012.fa
# 58592
grep -c "^>" concatenated_Clade_A_genomes.fasta
# 223223 # Concatenated transcript count equals sum of individual references.
```

Concatenate Algal genomes into Master Algal Genome
```
cat concatenated_Clade_A_genomes.fasta \
	concatenated_Clade_B_genomes.fasta \
	concatenated_Clade_C_genomes.fasta \
	concatenated_Clade_D_genomes_unique.fasta \
	> concatenated_Algal_genomes.fasta
```

Make Master Algal Symbiont Database

```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb \
    -in concatenated_Algal_genomes.fasta \
    -parse_seqids -dbtype nucl -out MasterAlgal_db
```

</details>

<details>

<summary>Make Bacteria Database</summary>

### Bacteria References:

I Followed the instructions in step 2 [here](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#:~:text=To%20use%20the%20download%20service,button%20to%20start%20the%20download)
Essentially, this downloaded 38,177 bacteria genomes (47.5 gb) that fit this criteria Bacteria, RefSeq, Complete Genomes. 
Ideally this can be used to blast our metatranscriptomes against to identify bacteria transcripts. 
A concise [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965038/) on the different available categories of bacteria genomes on NCBI: 
Not used: [Overview of DIAMOND and MEGAN](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.59)

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
```
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
Concatenate all the bacteria genomes into one multi-fasta. then use makeblastdb to make that multi-fasta a reference database.
```
scp ./genome_assemblies_genome_fasta.tar nicholas.macknight@holocron:/home/cns.local/nicholas.macknight/references/bacteria_reference/
```
> Runtime: 40 min. 
Decompress and concatenate with zcat
> /home/cns.local/nicholas.macknight/references/bacteria_reference
```
zcat *.fna.gz > concatenated_bacteria_genomes.fna
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/makeblastdb -in concatenated_bacteria_genomes.fna -parse_seqids -dbtype nucl -out MasterBacteria_db

### Acer
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/BacteriaOnly.Trinity.txt -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' BacteriaOnly.Trinity.txt > Acer_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Acer_Bacteria_contigs_percent_95.txt > Acer_Bacteria_contigs_percent_95_bp_150.txt

### Mcav
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/BacteriaOnly.Trinity.txt -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' LongestIsoform.BacteriaOnly.Trinity.txt > Mcav_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Mcav_Bacteria_contigs_percent_95.txt > Mcav_Bacteria_contigs_percent_95_bp_150.txt

### Ofav
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/trinity_out_dir.AllOfavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/BacteriaOnly.Trinity.txt -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' Ofav_trinity_output.LongestIsoform.BacteriaOnly.Trinity.txt > Ofav_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Ofav_Bacteria_contigs_percent_95.txt > Ofav_Bacteria_contigs_percent_95_bp_150.txt

### Past
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/BacteriaOnly.Trinity.txt -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' LongestIsoform.BacteriaOnly.Trinity.txt > Past_Bacteria_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Past_Bacteria_contigs_percent_95.txt > Past_Bacteria_contigs_percent_95_bp_150.txt

```

Now that we have a text file of the transcript names that pass our threshold for their alignment with the holobiont reference (example: Past_Bacteria_contigs_percent_95_bp_150.txt), we will extract the sequences from the metatranscriptomes (example: trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta) by first indexing the metatranscriptome with cdbfasta and then producing the .fa files (example: Past_Bacteria_only_transcriptome.fa).

:Microbe: Acer Bacteria
```
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta

cat Acer_Bacteria_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Acer_Bacteria_only_transcriptome.fa
```
> Run time: <1 min.

:microbe: Mcav Bacteria
```
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta

cat Mcav_Bacteria_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Mcav_Bacteria_only_transcriptome.fa
```

:microbe: Ofav Bacteria
```
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllOfavSamples_Lane1-8.LongestIsoform.Trinity.fasta

cat Ofav_Bacteria_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllOfavSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Ofav_Bacteria_only_transcriptome.fa
```

:microbe: Past Bacteria
```
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta

cat Past_Bacteria_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Past_Bacteria_only_transcriptome.fa
```

To count the number of transcripts in a fasta file
```
grep -c ">" filename.fa

# Example:
grep -c ">" Past_Bacteria_only_transcriptome.fa
```
</details>

<details>

<summary>Blastn - Assign metatranscriptome sequences to Host, Algal Symbiont, or Bacteria</summary>

**Note to self, consider removing host and bacteria blastn portions of this Blastn section as being redundant to previous code which is already doing the blastn. What is below may be a workshopped older version. The algal aspect of this section is up to date though.** 

🪸 Acer
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastx -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/trinity_out_dir_AllAcerSamples_Lane1-8/trinity_out_dir.LongestIsoform.CoralOnly.Trinity.txt
```
:microbe: Acer
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.LongestIsoform.BacteriaOnly.Trinity.txt -num_threads 20
```

🪸 Past 
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20
```

:microbe: Past
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/LongestIsoform.BacteriaOnly.Trinity.txt -num_threads 20
```

🪸 Ofav
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.CoralOnly.Trinity.txt -num_threads 20
```

:microbe: Ofav
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.BacteriaOnly.Trinity.txt -num_threads 20
```

🪸 Mcav
**blastn - Coral Host**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/Coral_Host/MasterCoral/MasterCoral_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.LongestIsoform.NewCoralOnly.Trinity.txt -num_threads 20
```

:microbe: Mcav
**blastn - Bacteria Only**
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta -db /home/cns.local/nicholas.macknight/references/bacteria_reference/MasterBacteria_db -outfmt "6 qseqid evalue pident length" -max_target_seqs 1 -out /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/BacteriaOnly.Trinity.txt -num_threads 20
```

### Blastn - Algal Symbionts

🌱 Acer - A
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_db/Clade_A_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_A_blastn_Acer_results.txt \
    -num_threads 20
```
🌱 Acer - B
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_B_db/Clade_B_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_B_blastn_Acer_results.txt \
    -num_threads 20
```
🌱 Acer - C
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_db/Clade_C_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_C_blastn_Acer_results.txt \
    -num_threads 20
```
🌱 Acer - D
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_db/Clade_D_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_D_blastn_Acer_results.txt \
    -num_threads 20
```
🌱 Mcav - A
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_db/Clade_A_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_A_blastn_Mcav_results.txt \
    -num_threads 20
```
🌱 Mcav - B
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_B_db/Clade_B_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_B_blastn_Mcav_results.txt \
    -num_threads 20
```
🌱 Mcav - C
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_db/Clade_C_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_C_blastn_Mcav_results.txt \
    -num_threads 20
```
🌱 Mcav - D
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_db/Clade_D_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_D_blastn_Mcav_results.txt \
    -num_threads 20
```
🌱 Ofav - A
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_db/Clade_A_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_A_blastn_Ofav_results.txt \
    -num_threads 20
```
🌱 Ofav - B
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_B_db/Clade_B_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_B_blastn_Ofav_results.txt \
    -num_threads 20
```
🌱 Ofav - C
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_db/Clade_C_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_C_blastn_Ofav_results.txt \
    -num_threads 20
```
🌱 Ofav - D
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_db/Clade_D_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_D_blastn_Ofav_results.txt \
    -num_threads 20

```
🌱 Past - A
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_db/Clade_A_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_A_blastn_Past_results.txt \
    -num_threads 20
```
🌱 Past - B
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_B_db/Clade_B_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_B_blastn_Past_results.txt \
    -num_threads 20
```
🌱 Past - C
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_db/Clade_C_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_C_blastn_Past_results.txt \
    -num_threads 20
```
🌱 Past - D
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_db/Clade_D_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out Clade_D_blastn_Past_results.txt \
    -num_threads 20
```

Algal Blastn to MastAlgalDB

🌱 Acer - Master Algal DB
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/MasterAlgal_db/MasterAlgal_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out MasterAlgal_blastn_Acer_results.txt \
    -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' MasterAlgal_blastn_Acer_results.txt > Acer_MasterAlgal_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Acer_MasterAlgal_contigs_percent_95.txt > Acer_MasterAlgal_contigs_percent_95_bp_150.txt

cat Acer_MasterAlgal_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Acer_MasterAlgal_only_transcriptome.fa

```
🌱 Mcav - Master Algal DB
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/MasterAlgal_db/MasterAlgal_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out MasterAlgal_blastn_Mcav_results.txt \
    -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' MasterAlgal_blastn_Mcav_results.txt > Mcav_MasterAlgal_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Mcav_MasterAlgal_contigs_percent_95.txt > Mcav_MasterAlgal_contigs_percent_95_bp_150.txt

cat Mcav_MasterAlgal_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Mcav_MasterAlgal_only_transcriptome.fa

```
🌱 Ofav - Master Algal DB
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_trinity_output.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/MasterAlgal_db/MasterAlgal_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out MasterAlgal_blastn_Ofav_results.txt \
    -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' MasterAlgal_blastn_Ofav_results.txt > Ofav_MasterAlgal_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Ofav_MasterAlgal_contigs_percent_95.txt > Ofav_MasterAlgal_contigs_percent_95_bp_150.txt

cat Ofav_MasterAlgal_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_trinity_output.LongestIsoform.Trinity.fasta.cidx > Ofav_MasterAlgal_only_transcriptome.fa

```
🌱 Past - Master Algal DB
```
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastn \
    -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta \
    -db /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/MasterAlgal_db/MasterAlgal_db \
    -outfmt "6 qseqid evalue pident length" \
    -max_target_seqs 1 \
    -out MasterAlgal_blastn_Past_results.txt \
    -num_threads 20

# Reads with less than 95% percent identity and shorter than 150 bp long are filtered out:
awk '{if ($3 > 95) print $1,$2,$4 }' MasterAlgal_blastn_Past_results.txt > Past_MasterAlgal_contigs_percent_95.txt
awk '{if ($3 > 150) print $1}' Past_MasterAlgal_contigs_percent_95.txt > Past_MasterAlgal_contigs_percent_95_bp_150.txt

cat Past_MasterAlgal_contigs_percent_95_bp_150.txt |  /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta.cidx > Past_MasterAlgal_only_transcriptome.fa

```
:seedling: :seedling: :seedling:

Threshold Filter contig blast quality - Loop "contig_evalue_filter.sh"

> What this script "contig_evalue_filter.sh" does is retain the contigs with percent identity > 95 AND contigs with length > 150 bp into a new .txt file. It then indexes the metranscriptome in preparation for the next script.
```
nano
```
copy and paiste (update file names to match yours and save as "contig_evalue_filter.sh"):
```
#!/bin/bash

# Define an array of clade names and coral hosts (modify as needed)
clades=("A" "B" "C" "D")
hosts=("Acer" "Mcav" "Ofav" "Past")

# Loop through each clade and host combination
for clade in "${clades[@]}"; do
    for host in "${hosts[@]}"; do
        # Construct the input and output filenames dynamically
        blastn_file="Clade_${clade}_blastn_${host}_results.txt"
        contigs_percent_95="Clade_${clade}_contigs_percent_95_${host}.txt"
        contigs_percent_95_bp_150="Clade_${clade}_contigs_percent_95_bp_150_${host}.txt"

        # Check if the blastn result file exists
        if [[ -f "$blastn_file" ]]; then
            echo "Processing $blastn_file..."

            # Extract contigs with percent identity > 95
            awk '{if ($3 > 95) print $1, $2, $4}' "$blastn_file" > "$contigs_percent_95"

            # Extract contigs with length > 150 bp
            awk '{if ($3 > 150) print $1}' "$contigs_percent_95" > "$contigs_percent_95_bp_150"
        else
            echo "Warning: $blastn_file not found. Skipping..."
        fi
    done
done

echo "Processing complete!"

# Create .cidx version with cdbfasta for contig extraction

/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta \
    /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta

/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta \
    trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta

/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta \
    trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta

/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta \
    /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_trinity_output.LongestIsoform.Trinity.fasta
```
:seedling: :seedling: :seedling:

Extract Contigs that pass Threshold - "extract_contigs.sh"

> This script takes the contig names from the previous script in the .txt file that passed the filtering thresholds and now this script extracts the sequences of those contig names into a new .fa file such as "Clade_A_Acer_algae_only_transcriptome.fa"
```
nano
```
Copy and paiste (update file names to match yours and save as "contig_evalue_filter.sh"):
```
#!/bin/bash

# Define arrays of clades and hosts
clades=("A" "B" "C" "D")
hosts=("Acer" "Mcav" "Past" "Ofav")

# Loop through clade and host combinations
for clade in "${clades[@]}"; do
    for host in "${hosts[@]}"; do

        # Set transcriptome paths and filenames based on the host
        if [[ "$host" == "Acer" ]]; then
            transcriptome_path="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer"
            transcriptome_file="trinity_out_dir.AllAcerSamples_Lane1-8.LongestIsoform.Trinity.fasta"
        elif [[ "$host" == "Mcav" ]]; then
            transcriptome_path="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav"
            transcriptome_file="trinity_out_dir.AllMcavSamples_Lane1-8.LongestIsoform.Trinity.fasta"
        elif [[ "$host" == "Past" ]]; then
            transcriptome_path="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past"
            transcriptome_file="trinity_out_dir.AllPastSamples_Lane1-8.LongestIsoform.Trinity.fasta"
        elif [[ "$host" == "Ofav" ]]; then
            transcriptome_path="/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output"
            transcriptome_file="Ofav_trinity_output.LongestIsoform.Trinity.fasta"
        fi

        # Construct full paths for the transcriptome and .cidx file
        transcriptome="${transcriptome_path}/${transcriptome_file}"
        cidx_file="${transcriptome}.cidx"

        # Construct filenames for filtered contigs and output fasta
        filtered_contigs="Clade_${clade}_contigs_percent_95_bp_150_${host}.txt"
        output_fasta="Clade_${clade}_${host}_algae_only_transcriptome.fa"

        # Check if both the filtered contigs and .cidx file exist
        if [[ -f "$filtered_contigs" && -f "$cidx_file" ]]; then
            echo "Extracting contigs for Clade $clade - $host..."

            # Extract contigs using cdbyank
            cat "$filtered_contigs" | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank "$cidx_file" > "$output_fasta"
        else
            echo "Warning: Missing $filtered_contigs or $cidx_file. Skipping..."
        fi
    done
done

echo "Contig extraction complete!"
```

:seedling: Concatenate Clade-Specific Transcripts Across Host
```
# Clade A
cat Clade_A_Acer_algae_only_transcriptome.fa \
    Clade_A_Mcav_algae_only_transcriptome.fa \
    Clade_A_Past_algae_only_transcriptome.fa \
    Clade_A_Ofav_algae_only_transcriptome.fa \
    > Master_Clade_A_transcriptome.fa

# Clade B
cat Clade_B_Acer_algae_only_transcriptome.fa \
    Clade_B_Mcav_algae_only_transcriptome.fa \
    Clade_B_Past_algae_only_transcriptome.fa \
    Clade_B_Ofav_algae_only_transcriptome.fa \
    > Master_Clade_B_transcriptome.fa

# Clade C
cat Clade_C_Acer_algae_only_transcriptome.fa \
    Clade_C_Mcav_algae_only_transcriptome.fa \
    Clade_C_Past_algae_only_transcriptome.fa \
    Clade_C_Ofav_algae_only_transcriptome.fa \
    > Master_Clade_C_transcriptome.fa

# Clade D
cat Clade_D_Acer_algae_only_transcriptome.fa \
    Clade_D_Mcav_algae_only_transcriptome.fa \
    Clade_D_Past_algae_only_transcriptome.fa \
    Clade_D_Ofav_algae_only_transcriptome.fa \
    > Master_Clade_D_transcriptome.fa
```

:seedling: Collapse to Longest Isoform: 
```
/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
    Master_Clade_A_transcriptome.fa > Master_Clade_A_longest_isoform.fa

/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
    Master_Clade_B_transcriptome.fa > Master_Clade_B_longest_isoform.fa

/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
    Master_Clade_C_transcriptome.fa > Master_Clade_C_longest_isoform.fa

/home/cns.local/nicholas.macknight/software/Trinity/trinityrnaseq-v2.15.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
    Master_Clade_D_transcriptome.fa > Master_Clade_D_longest_isoform.fa
```
</details>
<details>
<summary>Inferring Algal Composition</summary>

:seedling: :seedling: :seedling:

>When starting this project, I intended to infer algal composition among samples by using bbsplit and using mapping percentage to algal references as a way to infer algal composition. After much troubleshooting and improvement with that approach I switched methods and inferred algal composition through blastn results by counting the number of contigs. Let me explain why this decision was made.


>Originally, the idea was to have the four most common algal genera as references in bbsplit. I chose algal references others had recently published with and the results were that each species had a mixed algal composition, with no symbiont genera representing nearly all reads. Biologically, we know these coral hosts are typically dominated (+90%) with a single symbiont genera, however they can symbiont switch, but the likely-hood that was occuring in all my corals was unlikely. So began the journey to scrutinize this initial approach to ensure that my results were biological and not a technical artificat introduced somewhere. The first thing I did was question the references. New algal references are created with improved sequencing and bioinformatics tools. So there is the quality of the references, but also what if a algal reference came from one of the same coral hosts I was using, in a more similar geographic location as my project's samples and other algal references were from distantly related coral hosts from the Pacific? I would expect biological explained lower mapping rates. I tested for this, by taking nearly all known algal references on NCBI and having them as references in bbsplit, where I allowed bbsplit to map ambiguous reads multiple times and the results are a spectrum of mapping percentages. This is due to the technical and biological biases I mentioned. Here is a google excel that summarizes the quality of many available algal references. Anecdotally I was looking at the N50 (the higher the better) and the number of genes (Algal symbionts have a striking 120,000 genes due to evolutionary duplication events (typical plants), this is compared to the expected 20,000 in eukaryotic coral host genomes). To tackle the possible technical and biological bias I mentioned, I decided to lean towards the highest quality algal references rather than those that simply had the highest mapping rate which was my original plan. Collectively this led to a noticeable improvement in inferred algal composition. However I was still not content and after much troubleshooting with bbsplit parameters and ways to organize the references I trialed out a new method because of this thought: Expression does not equal abundance, a low abundance organism may appear overly abundant if its highly expressive and we still may be experiencing some technical bias factors I mentioned. So I decided to try to infer algal composition based on the number of unique contigs to remove the influence of expression over inflating inferred composition. The idea being that originally in bbsplit, if there were multiple reads of a contig, then that increases the mapping rate even though its the same contig. In this new approach once a unique contig has been identified, additional expression of that contig has no effect. The methods for this approach are written out in detail in this pipeline, but an overview is that I would blast the metranscriptomes to the algal references producing a .txt file of the transcript names that mapped. Through cdbyank, those transcript names were identified in the metatranscriptome, so that the sequences of those transcript names were extracted into the .fa files such as "Clade_A_Acer_algae_only_transcriptome.fa" then I could simply grep -c "^>" Clade_A_Acer_algae_only_transcriptome.fa to get a transcript count (contig count). I could see how many transcripts wer in Clade_B_Acer... Clade_C_Acer... and Clade_D_Acer... I would note the contig count and here is a summary of that:

<img width="1148" alt="image" src="https://github.com/user-attachments/assets/220b7ae0-921e-4391-81bf-57a2c73c56e8" />

|Species| BBsplit - Historical References  | BBsplit - Highest Quality References | Blastn - Number of Contigs |
| ------------- | ------------- | ------------- | ------------- |
| Acer | 38% - A  | 95% - A  | 99% - A |
| Mcav | 50% - B  | 68% - C  | 95% - C|
| Ofav | 70% - D  | 85% - D  | 99% - D|
| Past | 40% - A  | 93% - A  | 98% - A|


In summary, expression does not equal abundance, as we already knew, but by addressing this, we now have methods to circumnavigate that and produce results that reflect the true composition of these algal symbionts.


</details>

<details>
<summary>Transdecoder - Convert nucleotides into predicted proteins for OrthoFinder</summary>

> Transdecoder first appeared as a straight-forward "convert nucleotides to protein prediction" intermediate step. I quickly realized a significant overlooked issue with this that I have addressed in my methods below.

> The Problem: Orthofinder relies on the predicted proteins from Transdecoder. To predict proteins, transdecoder identifies Open Reading Frames "ORFs". Open reading frames begin with a start codon AUG and end with a stop codon and are the nucleotides in between which are potentially protein-coding regions. The problem is that there can be multiple, possibly hundreds of ORFs per gene. So, Orthofinder then using these multiple ORFs per gene to identify orthogroups, single copy orthologs, etc. These multiple ORF predicted orthogroups become a problem when its time to merge that information with count data produced from Salmon. Salmon uses the reference transcriptome that we used from BBsplit to count the frequency of those genes. So these are gene-level counts, NOT multiple ORF within a gene-level counts. So what do you do when you have lets say 750 ORFs, with unique functions, molecular evolution and need to assign a count at the gene-level? The answer is it doesnt work cleanly and creates a technical issue. There are various ways people address this and I played with most of them. You can pick the first ORF, the longest ORF, the ORF with the most homology. When I picked the first ORF, the problem is you have to just delete the other ORFs and that feels so messy and arbitrary. When I picked the longest ORF, I observed very strong species bias. This is probably because across mutliple species, some may evolutionarily have duplication events in their genes that will make them longer, or their references are higher quality leading to better alignment because its a better studied speices or their genome happened to be produced more recently with better equipment, or something about that species does better during the benchwork stage and extracts better than the others. When I was playing with the longest ORF approach I focused on the first function of transdecoder: transdecoder.LonOrfs, which states it selects the longest ORF. But when you read the details of this function, it does not actually retain ONLY the longest ORFs and has exceptions that will allow it keep multiple ORFs. I wrote a python script to force it to retain the longest ORF by identifying the transcript/gene names that were the same, and then keeping the one that was the longest length. These details can be found in the transcript titles. Once I did this it gave me the species bias I described and was not pursued any further in the analysis. I do not include the script to force the retention of longest ORF in this pipeline because I do not want to confuse any user that it is an approach I applied. The third approach, I ended up sticking with is a ORF homology-based approach. During the Transdecoder.Predict function, you can specify "single_best_only". This function will select only one representative ORF per transcript based on homology. During the prediction process, TransDecoder evaluates ORFs based on coding potential, length, and homology to known protein domains or sequences.
ORFs are ranked by these scores, with higher scores indicating a stronger likelihood of representing true protein-coding sequences as the "single-best" ORF. I found this homolgy based approach to represent the biological nature of the data while the others introduced a technical bias. 

> When not to select one ORF: Selecting one ORF is ideal for those who are interested in gene or transcript-level analysis and do not care about isoform or ORF level analysis. Much of multi-species comparisons is a tradeoff between complexity and comparability. If you want to do gene-level comparisons, there is likely to be more conservation of those genes and they are less species-specific and therefore more comparable. If you want to do ORF-level comparisons, you will retain the complexity of each species, but that complexity may be species-specific and less evolutionarity conserved. Because we are comparing four species in this study that expand potentially 100+ million years of evolutionary divergence, I am leaning into comparability. But if I was compariing Orbicella faveolata to Orbicella franksi or even Orbicella faveolata to Montastraea cavernosa (15 million year divergence), I may lean into retaining that complexity to explore the molecular evolution of these relatively closely related species. 

### *Retaining Non-Protein-Coding Regions*

> BONUS. From this deep dive into protein-prediction, I realized the default transdecoder only predicts protein-coding regions and excludes non-protein-coding regions. Why should we care about non-protein-coding regions? Non-protein-coding regions are highly valuable in multi-species RNA expression comparison studies, particularly in the context of disease response, because they often *regulate* gene expression, mediate cellular signaling, and adaptively respond/are hijacked to environmental or by pathological stimuli. Because our questions evolve around what is SCTLD and why are some species resistant, these non-protein-coding regions are worth retaining and this pipeline integrates how to do this.
To briefly explain, in adition to protein-coding protein prediction with transdecoder, we also perform pfam and a blastp search. This gives us three "datasets" of semi-overlapping results that we feed into TransDecoder.Predict and have it retain the "single_best_only" per transcript. When you just retain protein-coding regions, you are inadvertenly discarding about ~60% of the transcriptome. So implementing pfam and blastp allows you to maximize your dataset. 

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

### :microbe: Bacteria

:microbe: Acer Bacteria 
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


:microbe: Mcav Bacteria
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

:microbe: Ofav Bacteria
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

:microbe: Past Bacteria
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

### Host 🪸

🪸 Acer Host
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

🪸 Mcav Host
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

🪸 Ofav Host
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

🪸 Past Host
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
### :seedling: Algae

:seedling: Acer Algae
```
/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Clade_A_Acer_algae_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Clade_A_Acer_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blastp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Clade_A_Acer_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs/Clade_A_Acer_algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Clade_A_Acer_algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Clade_A_Acer_algae_only_transcriptome.fa.transdecoder.pep Clade_A_Acer_algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_A_Acer_algae_only_transcriptome_transdecoder.fa -o Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly.fa
```

:seedling: Mcav Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Clade_C_Mcav_algae_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Clade_C_Mcav_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Clade_C_Mcav_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs/Clade_C_Mcav_algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Clade_C_Mcav_algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Clade_C_Mcav_algae_only_transcriptome.fa.transdecoder.pep Clade_C_Mcav_algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_C_Mcav_algae_only_transcriptome_transdecoder.fa -o Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly.fa
```

:seedling: Ofav Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Clade_D_Ofav_algae_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs/Clade_D_Ofav_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs/Clade_D_Ofav_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs/Clade_D_Ofav_algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Clade_D_Ofav_algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Clade_D_Ofav_algae_only_transcriptome.fa.transdecoder.pep Clade_D_Ofav_algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_D_Ofav_algae_only_transcriptome_transdecoder.fa -o Clade_D_Ofav_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Clade_D_Ofav_reference_proteome_AllORF_SingleBestOnly.fa
```

:seedling: Past Algae
```
mkdir Algae_transdecoder_AllORFs
cd Algae_transdecoder_AllORFs

# TransDecoder.LongOrfs (Protein-Coding ORF Prediciton)
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ../Clade_A_Past_algae_only_transcriptome.fa

# Pfam
/home/cns.local/nicholas.macknight/software/hmmer-3.4/src/hmmsearch --cpu 50 -E 1e-10 --domtblout pfam.domtblout /home/cns.local/nicholas.macknight/software/Pfam-A.hmm.gz /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Clade_A_Past_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep

# Blasp
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Clade_A_Past_algae_only_transcriptome.fa.transdecoder_dir/longest_orfs.pep  \
    -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 50 > /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs/Clade_A_Past_algae_only_transcriptome.fa.transdecoder_dir/blastp.outfmt6

# Transdecoder.Predict 
/home/cns.local/nicholas.macknight/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ../Clade_A_Past_algae_only_transcriptome.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only

# CD-HIT
mv Clade_A_Past_algae_only_transcriptome.fa.transdecoder.pep Clade_A_Past_algae_only_transcriptome_transdecoder.fa
/home/cns.local/nicholas.macknight/software/cd-hit-v4.8.1-2019-0228/cd-hit -i Clade_A_Past_algae_only_transcriptome_transdecoder.fa -o Clade_A_Past_reference_proteome_AllORF_SingleBestOnly.fa

# Verify Single ORF
./check_single_orf.py Clade_A_Past_reference_proteome_AllORF_SingleBestOnly.fa
```
</details>
<details>
<summary>OrthoFinder - Identify homologous sequences among species</summary>

# OrthoFinder

Move all the reference_proteome.fa into a new folder specific for each holobiont compartment (Host, Algae, Bacteria, so three folders total)
Example of moving all reference_proteome.fa files into their respective newly created folder:

🪸 Host
```
# Acer
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF
scp /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Host_transdecoder_AllORFs/Acer_Host_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/

# Mcav
scp /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Host_transdecoder_AllORFs/Mcav_coral_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
# Mcav Rename
mv Mcav_coral_reference_proteome_AllORF_SingleBestOnly.fa Mcav_Host_reference_proteome_AllORF_SingleBestOnly.fa

#Ofav
scp /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Host_transdecoder_AllORFs/Ofav_Host_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/

# Past
scp /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Host_transdecoder_AllORFs/Past_Host_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
```

:seedling: Algae
```
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF

#Acer
cd /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Algae_transdecoder_AllORFs
scp Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/

#Mcav
cd /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Algae_transdecoder_AllORFs
scp Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/

#Ofav
cd /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Algae_transdecoder_AllORFs
scp Clade_A_Ofav_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/

#Past
cd /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Algae_transdecoder_AllORFs
scp Clade_A_Past_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/
```

:microbe: Bacteria
```
mkdir /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF

# Acer
# In performing this I used the mv command when I should have used scp. Keeping mv for records and making this note.
mv /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Bacteria_transdecoder_AllORFs/Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/

# Mcav
# In performing this I used the mv command when I should have used scp. Keeping mv for records and making this note.
mv /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Bacteria_transdecoder_AllORFs/Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/

# Ofav
# In performing this I used the mv command when I should have used scp. Keeping mv for records and making this note.
mv /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Bacteria_transdecoder_AllORFs/Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/

# Past
# In performing this I used the mv command when I should have used scp. Keeping mv for records and making this note.
mv /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Bacteria_transdecoder_AllORFs/Past_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/
```

Then within each folder as the active directory perform the Orthofinder Command:
> This command performs OrthoFinder, comparing the predicted proteins among the reference_proteome.fa files
```
# 🪸 Host
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50

# :seedling: Algae
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50

# 🦠 Bacteria
cd /home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/
python /home/cns.local/nicholas.macknight/software/OrthoFinder_source/orthofinder.py -f . -t 50
```

</details>

<details>
<summary>Annotating OrthoGroups</summary>
# Annotating OrthoGroups

Orthogroups contain multiple transcripts. It is common for these transcripts to have unique gene annotations. So how do researchers decide how to annotate the Orthogroup when there are unique transcript-level annotations to choose from? There are several approaches. I will begin with the approach I am applying (Most Common Annotation ID Approach) and you can choose to read about the others if you'd like. 

**Most Common Annotation ID Approach**: Simply put I choose the annotation that is the most common within that Orthogroup. So if there are ten transcripts within an Orthogroup and seven of them are the same, that that becomes the majority annotation and becomes the representative annotation for that Orthogroup. I chose this approach because to me it seems the most fair and does reduce reference bias in annotation. 

**Highest Quality Transcriptome Approach**: If you are performing multi species comparison transcriptomics, you will have multiple transcriptomes and references used to annotate those transcriptomes. This approach selects the species with the highest quality transcriptome, annotates those transcripts, and uses those transcript annotations to annotate the Orthogroup. This is a logical idea as the sentiment is your transcript annoations will be highest possible accuracy. The caveat is your annotations could be considered organism-biased.

**Quality Score Approach**: Each transcript annotation has a bit score (strength of alignment) and E-Value (which considered bit score plus signficance of alignment). Acceptable annotation cutoffs are field specific but for coral it is any annotation with an E-Value less than E-5 is acceptable. This annotation approach will select the E-Value with the lowest value as the representative annotation for that Orthgroup. This is also a sensible approach as you are picking the annotation that is the strongest. The caveat here is if we think of e values like p values, whats effectively the difference between a p value of lets say 0.00001 and 0.0000001. The mathematic scale of magnitude is of unique value but my interpretation and value that I assign it in my discussion is going to be the same. 

> Annotation organization occurs through command-line to assign the annotation name to each coral dataset. Then the preffered approach listed above is applied in R.


### :microbe: Annotating Bacteria Orthogroups 

First we need to create a file of the orthogroup transcripts. "_bac_Orthogroup_Transcripts.txt" from the Orthofinder results.
Move these files from the server to your local computer.
```
# On my local computer (not in the server):
cd /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/OrthoFinder/Results_Dec21/Orthogroups/Orthogroups.tsv ./
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/OrthoFinder/Results_Dec21/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/OrthoFinder/Results_Dec21/Orthogroups/Comparative_Genomics_Statistics ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/OrthoFinder/Results_Dec21/Orthogroups/Species_Tree ./
```
🦠 Organize Orthogroups.tsv to make "_bac_Orthogroup_Transcripts.txt" in R

> This can be done in excel like so:
> 1. Open Orthogroups.tsv in Excel
> 2. Click *Find & Select*
> 3. Click *Go to Special*
> 4. Choose *Blanks*
> 5. Click OK and then all the *blank rows/cells will be highlighted*
> 6. Choose the *Delete under Cells* section on the Home Tab
> 7. Click *Delete Sheet Rows*
> 8. Save as "Bacteria_shared_orthogroups_SingleBestORF.csv"
> Steps 1-8 were repeated 4 times so that all blank rows were eventually removed. I think
the amount of data excel needed to process required these steps to be
repeated.
You can ensure you have the accurate final number of
orthogroups by comparing your row count to the "Statistics_Overall" file
in Comparative_Genomics_Statistics specifically the value in row "Number of orthogroups with all species present

In R:
```
```{r Bacteria Orthogroups}

Orthogroups <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/Bacteria_shared_orthogroups_SingleBestORF.csv")

# Acer Orthogroups
Acer_orthogroups <- Orthogroups[,c(1,2)]
Acer_orthogroups <- Acer_orthogroups %>% separate_rows(Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Acer_orthogroups <- Acer_orthogroups[,c(2,1)]
colnames(Acer_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Acer_bac_Orthogroup_Transcripts <- Acer_orthogroups[c(1)]
names(Acer_bac_Orthogroup_Transcripts)
write.table(Acer_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Acer_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Acer_orthogroups$Transcript <- gsub("\\..*","",Acer_orthogroups$Transcript)
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Mcav Orthogroups
Mcav_orthogroups <- Orthogroups[,c(1,3)]
Mcav_orthogroups <- Mcav_orthogroups %>% separate_rows(Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Mcav_orthogroups <- Mcav_orthogroups[,c(2,1)]
colnames(Mcav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Mcav_bac_Orthogroup_Transcripts <- Mcav_orthogroups[c(1)]
names(Mcav_bac_Orthogroup_Transcripts)
write.table(Mcav_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Mcav_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Mcav_orthogroups$Transcript <- gsub("\\..*","",Mcav_orthogroups$Transcript)
Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/McavBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Ofav Orthogroups
Ofav_orthogroups <- Orthogroups[,c(1,4)]
Ofav_orthogroups <- Ofav_orthogroups %>% separate_rows(Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Ofav_orthogroups <- Ofav_orthogroups[,c(2,1)]
colnames(Ofav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Ofav_bac_Orthogroup_Transcripts <- Ofav_orthogroups[c(1)]
names(Ofav_bac_Orthogroup_Transcripts)
write.table(Ofav_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Ofav_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Ofav_orthogroups$Transcript <- gsub("\\..*","",Ofav_orthogroups$Transcript)
Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene.csv",row.names = FALSE)

# Past Orthogroups
Past_orthogroups <- Orthogroups[,c(1,5)]
Past_orthogroups <- Past_orthogroups %>% separate_rows(Past_Bacteria_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Past_orthogroups <- Past_orthogroups[,c(2,1)]
colnames(Past_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Past_bac_Orthogroup_Transcripts <- Past_orthogroups[c(1)]
names(Past_bac_Orthogroup_Transcripts)
write.table(Past_bac_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Bacteria/Past_bac_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Past_orthogroups$Transcript <- gsub("\\..*","",Past_orthogroups$Transcript)
Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/PastBac_orthogroup_tx2gene.csv",row.names = FALSE)

```

🦠 Annotating Transcripts
```
mkdir Annotating_Orthogroups
# Make an index of the reference proteome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Past_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa

# Extract the sequences of orthogroup transcripts from index: 
cat Acer_bac_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Acer_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Acer_Bac_orthologs.fa
cat Mcav_bac_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Mcav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Mcav_Bac_orthologs.fa
cat Ofav_bac_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Ofav_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Ofav_Bac_orthologs.fa
cat Past_bac_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Past_Bacteria_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Past_Bac_orthologs.fa

# Annotate extracted sequences with BLASTp:
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Acer_Bac_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Acer_Bac_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Mcav_Bac_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Mcav_Bac_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Ofav_Bac_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Ofav_Bac_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Past_Bac_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Past_Bac_orthologs_annotated.txt -num_threads 60


### Retain Only those with an e value of e^-5 or less
awk '{if ($3 < 1e-05) print $1,$2,$3}' Acer_Bac_orthologs_annotated.txt > Acer_Bac_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Mcav_Bac_orthologs_annotated.txt > Mcav_Bac_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Ofav_Bac_orthologs_annotated.txt > Ofav_Bac_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Past_Bac_orthologs_annotated.txt > Past_Bac_orthologs_annotated_e-5.txt

```
In R:
🦠 Read in Annotated Transcripts
```
# In terminal: 
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Bacteria/SingleBestORF/*_Bac_orthologs_annotated_e-5.txt /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/R/Bacteria/
  
# Load required libraries
library(dplyr)
library(readr)

# List of input file prefixes
file_prefixes <- c("Acer", "Mcav", "Ofav", "Past")

# Function to process a single file
process_file <- function(prefix) {
  # Construct file paths
  input_file <- paste0(prefix, "_Bac_orthologs_annotated_e-5.txt")
  output_file <- paste0(prefix, "_Bac_orthologs_annotated_e-5_formatted.txt")
  
  # Check if the input file exists
  if (!file.exists(input_file)) {
    message(paste("File not found:", input_file))
    return(NULL)
  }
  
  # Load the data
  data <- read.table(input_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
  
  # Format the data
  formatted_data <- data %>%
    mutate(Entry = sub(".*\\|(.*)\\|.*", "\\1", V1)) %>%
    select(Entry, Transcript = V2, Evalue = V3)
  
  # Write the output file
  write.table(formatted_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Return the formatted data for verification
  return(formatted_data)
}

# Process each file and store results in a list
formatted_results <- lapply(file_prefixes, process_file)

# Optionally, print the first few rows of each processed file
names(formatted_results) <- file_prefixes
lapply(formatted_results, function(x) {
  if (!is.null(x)) head(x) else NULL
})

Acer_annot_transcripts <- formatted_results$Acer
Acer_annot_transcripts$Transcript <- gsub("\\..*","",Acer_annot_transcripts$Transcript) # Remove ORF Identifier

Mcav_annot_transcripts <- formatted_results$Mcav
Mcav_annot_transcripts$Transcript <- gsub("\\..*","",Mcav_annot_transcripts$Transcript)

Ofav_annot_transcripts <- formatted_results$Ofav
Ofav_annot_transcripts$Transcript <- gsub("\\..*","",Ofav_annot_transcripts$Transcript)

Past_annot_transcripts <- formatted_results$Past
Past_annot_transcripts$Transcript <- gsub("\\..*","",Past_annot_transcripts$Transcript)


```
🦠 Annotate Orthogroups
```
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/AcerBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/McavBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/OfavBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Bacteria/SingleBestORF/PastBac_orthogroup_tx2gene_annot.csv",row.names = FALSE)


```

### :seedling: Annotating Algae Orthogroups 
First we need to create a file of the orthogroup transcripts. "_algae_Orthogroup_Transcripts.txt" from the Orthofinder results.
Move these files from the server to your local computer.
```
# On my local computer (not in the server):
cd /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/Algae/
mkdir Results_Jan13
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/OrthoFinder/Results_Jan13/Orthogroups/Orthogroups.tsv ./
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/OrthoFinder/Results_Jan13/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/OrthoFinder/Results_Jan13/Comparative_Genomics_Statistics ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Algae/SingleBestORF/OrthoFinder/Results_Jan13/Species_Tree ./
```
:seedling: Organize Orthogroups.tsv to make "_algae_Orthogroup_Transcripts.txt" in R

> This can be done in excel like so:
> 1. Open Orthogroups.tsv in Excel
> 2. Click *Find & Select*
> 3. Click *Go to Special*
> 4. Choose *Blanks*
> 5. Click OK and then all the *blank rows/cells will be highlighted*
> 6. Choose the *Delete under Cells* section on the Home Tab
> 7. Click *Delete Sheet Rows*
> 8. Save as "algae_shared_orthogroups_SingleBestORF.csv"
> Steps 1-8 were repeated 4 times so that all blank rows were eventually removed. I think
the amount of data excel needed to process required these steps to be
repeated.
You can ensure you have the accurate final number of
orthogroups by comparing your row count to the "Statistics_Overall" file
in Comparative_Genomics_Statistics specifically the value in row "Number of orthogroups with all species present

In R:
```
```{r algae Orthogroups}

Orthogroups <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Algae/Results_Jan13/algae_shared_orthogroups.csv")

# Acer Orthogroups
Acer_orthogroups <- Orthogroups[,c(1,2)]
Acer_orthogroups <- Acer_orthogroups %>% separate_rows(Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Acer_orthogroups <- Acer_orthogroups[,c(2,1)]
colnames(Acer_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Acer_algae_Orthogroup_Transcripts <- Acer_orthogroups[c(1)]
names(Acer_algae_Orthogroup_Transcripts)
write.table(Acer_algae_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/algae/Acer_algae_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Acer_orthogroups$Transcript <- gsub("\\..*","",Acer_orthogroups$Transcript)
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Clade_A_Acer_orthogroup_tx2gene_SingleBestOnly.csv",row.names = FALSE)

# Mcav Orthogroups
Mcav_orthogroups <- Orthogroups[,c(1,3)]
Mcav_orthogroups <- Mcav_orthogroups %>% separate_rows(Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Mcav_orthogroups <- Mcav_orthogroups[,c(2,1)]
colnames(Mcav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Mcav_algae_Orthogroup_Transcripts <- Mcav_orthogroups[c(1)]
names(Mcav_algae_Orthogroup_Transcripts)
write.table(Mcav_algae_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/algae/Mcav_algae_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Mcav_orthogroups$Transcript <- gsub("\\..*","",Mcav_orthogroups$Transcript)
Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Clade_C_Mcav_orthogroup_tx2gene_SingleBestOnly.csv",row.names = FALSE)

# Ofav Orthogroups
Ofav_orthogroups <- Orthogroups[,c(1,4)]
Ofav_orthogroups <- Ofav_orthogroups %>% separate_rows(Clade_D_Ofav_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Ofav_orthogroups <- Ofav_orthogroups[,c(2,1)]
colnames(Ofav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Ofav_algae_Orthogroup_Transcripts <- Ofav_orthogroups[c(1)]
names(Ofav_algae_Orthogroup_Transcripts)
write.table(Ofav_algae_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/algae/Ofav_algae_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Ofav_orthogroups$Transcript <- gsub("\\..*","",Ofav_orthogroups$Transcript)
Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Clade_D_Ofav_orthogroup_tx2gene_SingleBestOnly.csv",row.names = FALSE)

# Past Orthogroups
Past_orthogroups <- Orthogroups[,c(1,5)]
Past_orthogroups <- Past_orthogroups %>% separate_rows(Clade_A_Past_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Past_orthogroups <- Past_orthogroups[,c(2,1)]
colnames(Past_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Past_algae_Orthogroup_Transcripts <- Past_orthogroups[c(1)]
names(Past_algae_Orthogroup_Transcripts)
write.table(Past_algae_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/algae/Past_algae_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Past_orthogroups$Transcript <- gsub("\\..*","",Past_orthogroups$Transcript)
Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Clade_A_Past_orthogroup_tx2gene_SingleBestOnly.csv",row.names = FALSE)

```
:seedling: Transfer Orthogroup Transcripts from local computer to server
```
# on local computer terminal
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Acer_algae_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/algae/Results_Jan13/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Mcav_algae_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/algae/Results_Jan13/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Ofav_algae_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/algae/Results_Jan13/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Past_algae_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/algae/Results_Jan13/


```

:seedling: Annotating Transcripts
```
mkdir Annotating_Orthogroups
# Make an index of the reference proteome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Clade_D_Ofav_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Clade_A_Past_reference_proteome_AllORF_SingleBestOnly.fa

# Extract the sequences of orthogroup transcripts from index: 
cat Acer_algae_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Clade_A_Acer_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Acer_algae_orthologs.fa
cat Mcav_algae_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Clade_C_Mcav_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Mcav_algae_orthologs.fa
cat Ofav_algae_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Clade_D_Ofav_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Ofav_algae_orthologs.fa
cat Past_algae_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Clade_A_Past_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Past_algae_orthologs.fa

# Annotate extracted sequences with BLASTp:
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Acer_algae_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Acer_algae_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Mcav_algae_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Mcav_algae_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Ofav_algae_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Ofav_algae_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Past_algae_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Past_algae_orthologs_annotated.txt -num_threads 60


### Retain Only those with an e value of e^-5 or less
awk '{if ($3 < 1e-05) print $1,$2,$3}' Acer_algae_orthologs_annotated.txt > Acer_algae_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Mcav_algae_orthologs_annotated.txt > Mcav_algae_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Ofav_algae_orthologs_annotated.txt > Ofav_algae_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Past_algae_orthologs_annotated.txt > Past_algae_orthologs_annotated_e-5.txt

```
In R:
:seedling: Read in Annotated Transcripts
```
# In terminal: 
scp nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/algae/Results_Jan13/*_algae_orthologs_annotated_e-5.txt /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/R/Algae/
  
# Load required libraries
library(dplyr)
library(readr)

# List of input file prefixes
file_prefixes <- c("Acer", "Mcav", "Ofav", "Past")

# Function to process a single file
process_file <- function(prefix) {
  # Construct file paths
  input_file <- paste0(prefix, "_algae_orthologs_annotated_e-5.txt")
  output_file <- paste0(prefix, "_algae_orthologs_annotated_e-5_formatted.txt")
  
  # Check if the input file exists
  if (!file.exists(input_file)) {
    message(paste("File not found:", input_file))
    return(NULL)
  }
  
  # Load the data
  data <- read.table(input_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
  
  # Format the data
  formatted_data <- data %>%
    mutate(Entry = sub(".*\\|(.*)\\|.*", "\\1", V1)) %>%
    select(Entry, Transcript = V2, Evalue = V3)
  
  # Write the output file
  write.table(formatted_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Return the formatted data for verification
  return(formatted_data)
}

# Process each file and store results in a list
formatted_results <- lapply(file_prefixes, process_file)

# Optionally, print the first few rows of each processed file
names(formatted_results) <- file_prefixes
lapply(formatted_results, function(x) {
  if (!is.null(x)) head(x) else NULL
})

Acer_annot_transcripts <- formatted_results$Acer
Acer_annot_transcripts$Transcript <- gsub("\\..*","",Acer_annot_transcripts$Transcript) # Remove ORF Identifier

Mcav_annot_transcripts <- formatted_results$Mcav
Mcav_annot_transcripts$Transcript <- gsub("\\..*","",Mcav_annot_transcripts$Transcript)

Ofav_annot_transcripts <- formatted_results$Ofav
Ofav_annot_transcripts$Transcript <- gsub("\\..*","",Ofav_annot_transcripts$Transcript)

Past_annot_transcripts <- formatted_results$Past
Past_annot_transcripts$Transcript <- gsub("\\..*","",Past_annot_transcripts$Transcript)


```
:seedling: Annotate Orthogroups
```
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Aceralgae_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Mcavalgae_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Ofavalgae_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/algae/Results_Jan13/Pastalgae_orthogroup_tx2gene_annot.csv",row.names = FALSE)


```

### 🪸 Annotating Host Orthogroups 
First we need to create a file of the orthogroup transcripts. "_host_Orthogroup_Transcripts.txt" from the Orthofinder results.
Move these files from the server to your local computer.
```
# On my local computer (not in the server):
cd /Users/nicholas.macknight/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/
mkdir Results_Jan13
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/OrthoFinder/Results_Jan10/Orthogroups/Orthogroups.tsv ./
scp nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/OrthoFinder/Results_Jan10/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/OrthoFinder/Results_Jan10/Comparative_Genomics_Statistics ./
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/OrthoFinder/Results_Jan10/Species_Tree ./
```


🪸 Organize Orthogroups.tsv to make "_host_Orthogroup_Transcripts.txt" in R

> This can be done in excel like so:
> 1. Open Orthogroups.tsv in Excel
> 2. Click *Find & Select*
> 3. Click *Go to Special*
> 4. Choose *Blanks*
> 5. Click OK and then all the *blank rows/cells will be highlighted*
> 6. Choose the *Delete under Cells* section on the Home Tab
> 7. Click *Delete Sheet Rows*
> 8. Save as "Host_shared_orthogroups_SingleBestORF.csv"
> Steps 1-8 were repeated 4 times so that all blank rows were eventually removed. I think
the amount of data excel needed to process required these steps to be
repeated.
You can ensure you have the accurate final number of
orthogroups by comparing your row count to the "Statistics_Overall" file
in Comparative_Genomics_Statistics specifically the value in row "Number of orthogroups with all species present

In R:
```
```{r Host Orthogroups}

Orthogroups <- read.csv("~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/host_shared_orthogroups_SingleBestORF.csv")

# Acer Orthogroups
Acer_orthogroups <- Orthogroups[,c(1,2)]
Acer_orthogroups <- Acer_orthogroups %>% separate_rows(Acer_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Acer_orthogroups <- Acer_orthogroups[,c(2,1)]
colnames(Acer_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Acer_host_Orthogroup_Transcripts <- Acer_orthogroups[c(1)]
names(Acer_host_Orthogroup_Transcripts)
write.table(Acer_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Host/Acer_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Acer_orthogroups$Transcript <- gsub("\\..*","",Acer_orthogroups$Transcript)
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Acerhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Mcav Orthogroups
Mcav_orthogroups <- Orthogroups[,c(1,3)]
Mcav_orthogroups <- Mcav_orthogroups %>% separate_rows(Mcav_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Mcav_orthogroups <- Mcav_orthogroups[,c(2,1)]
colnames(Mcav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Mcav_host_Orthogroup_Transcripts <- Mcav_orthogroups[c(1)]
names(Mcav_host_Orthogroup_Transcripts)
write.table(Mcav_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Host/Mcav_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Mcav_orthogroups$Transcript <- gsub("\\..*","",Mcav_orthogroups$Transcript)
Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Mcavhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Ofav Orthogroups
Ofav_orthogroups <- Orthogroups[,c(1,4)]
Ofav_orthogroups <- Ofav_orthogroups %>% separate_rows(Ofav_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Ofav_orthogroups <- Ofav_orthogroups[,c(2,1)]
colnames(Ofav_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Ofav_host_Orthogroup_Transcripts <- Ofav_orthogroups[c(1)]
names(Ofav_host_Orthogroup_Transcripts)
write.table(Ofav_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Host/Ofav_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Ofav_orthogroups$Transcript <- gsub("\\..*","",Ofav_orthogroups$Transcript)
Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Ofavhost_orthogroup_tx2gene.csv",row.names = FALSE)

# Past Orthogroups
Past_orthogroups <- Orthogroups[,c(1,5)]
Past_orthogroups <- Past_orthogroups %>% separate_rows(Past_Host_reference_proteome_AllORF_SingleBestOnly, sep = ", ")
Past_orthogroups <- Past_orthogroups[,c(2,1)]
colnames(Past_orthogroups)[1] <- "Transcript"

### Pull out Orthogroup transcripts for annotation:
Past_host_Orthogroup_Transcripts <- Past_orthogroups[c(1)]
names(Past_host_Orthogroup_Transcripts)
write.table(Past_host_Orthogroup_Transcripts, file="~/Desktop/Microbial Metatranscriptomics/R/Host/Past_host_Orthogroup_Transcripts.txt",quote = FALSE,row.names = FALSE)

### Format Orthogroup Transcripts For R: 
Past_orthogroups$Transcript <- gsub("\\..*","",Past_orthogroups$Transcript)
Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Pasthost_orthogroup_tx2gene.csv",row.names = FALSE)

```
🪸 Transfer Orthogroup Transcripts from local computer to server
```
# on local computer terminal
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Acer_host_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Mcav_host_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Ofav_host_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/
scp ~/Desktop/Microbial\ Metatranscriptomics/OrthoFinder/Host/Results_Jan10/Past_host_Orthogroup_Transcripts.txt nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/


```

🪸 Annotating Transcripts
```
mkdir Annotating_Orthogroups
# Make an index of the reference proteome
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Acer_Host_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Mcav_Host_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Ofav_Host_reference_proteome_AllORF_SingleBestOnly.fa
/home/cns.local/nicholas.macknight/software/cdbfasta/cdbfasta Past_Host_reference_proteome_AllORF_SingleBestOnly.fa

# Extract the sequences of orthogroup transcripts from index: 
cat Acer_host_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Acer_Host_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Acer_host_orthologs.fa
cat Mcav_host_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Mcav_Host_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Mcav_host_orthologs.fa
cat Ofav_host_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Ofav_Host_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Ofav_host_orthologs.fa
cat Past_host_Orthogroup_Transcripts.txt | /home/cns.local/nicholas.macknight/software/cdbfasta/cdbyank Past_Host_reference_proteome_AllORF_SingleBestOnly.fa.cidx > Past_host_orthologs.fa

# Annotate extracted sequences with BLASTp:
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Acer_host_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Acer_host_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Mcav_host_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Mcav_host_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Ofav_host_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Ofav_host_orthologs_annotated.txt -num_threads 60
/home/cns.local/nicholas.macknight/software/ncbi-blast-2.15.0+/bin/blastp -query Past_host_orthologs.fa -db /home/cns.local/nicholas.macknight/references/uniprot/uniprot_db -outfmt "6 sseqid qseqid evalue" -max_target_seqs 1 -out Past_host_orthologs_annotated.txt -num_threads 60


### Retain Only those with an e value of e^-5 or less
awk '{if ($3 < 1e-05) print $1,$2,$3}' Acer_host_orthologs_annotated.txt > Acer_host_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Mcav_host_orthologs_annotated.txt > Mcav_host_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Ofav_host_orthologs_annotated.txt > Ofav_host_orthologs_annotated_e-5.txt
awk '{if ($3 < 1e-05) print $1,$2,$3}' Past_host_orthologs_annotated.txt > Past_host_orthologs_annotated_e-5.txt

```
In R:
🪸 Read in Annotated Transcripts
```
# In terminal: 
scp nicholas.macknight@holocron:../../../home/cns.local/nicholas.macknight/SCTLDRNA/Orthofinder/Host/SingleBestORF/*_host_orthologs_annotated_e-5.txt /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/R/Coral\ Host/
  
# Load required libraries
library(dplyr)
library(readr)

# List of input file prefixes
file_prefixes <- c("Acer", "Mcav", "Ofav", "Past")

# Function to process a single file
process_file <- function(prefix) {
  # Construct file paths
  input_file <- paste0(prefix, "_host_orthologs_annotated_e-5.txt")
  output_file <- paste0(prefix, "_host_orthologs_annotated_e-5_formatted.txt")
  
  # Check if the input file exists
  if (!file.exists(input_file)) {
    message(paste("File not found:", input_file))
    return(NULL)
  }
  
  # Load the data
  data <- read.table(input_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
  
  # Format the data
  formatted_data <- data %>%
    mutate(Entry = sub(".*\\|(.*)\\|.*", "\\1", V1)) %>%
    select(Entry, Transcript = V2, Evalue = V3)
  
  # Write the output file
  write.table(formatted_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Return the formatted data for verification
  return(formatted_data)
}

# Process each file and store results in a list
formatted_results <- lapply(file_prefixes, process_file)

# Optionally, print the first few rows of each processed file
names(formatted_results) <- file_prefixes
lapply(formatted_results, function(x) {
  if (!is.null(x)) head(x) else NULL
})

Acer_annot_transcripts <- formatted_results$Acer
Acer_annot_transcripts$Transcript <- gsub("\\..*","",Acer_annot_transcripts$Transcript) # Remove ORF Identifier

Mcav_annot_transcripts <- formatted_results$Mcav
Mcav_annot_transcripts$Transcript <- gsub("\\..*","",Mcav_annot_transcripts$Transcript)

Ofav_annot_transcripts <- formatted_results$Ofav
Ofav_annot_transcripts$Transcript <- gsub("\\..*","",Ofav_annot_transcripts$Transcript)

Past_annot_transcripts <- formatted_results$Past
Past_annot_transcripts$Transcript <- gsub("\\..*","",Past_annot_transcripts$Transcript)


```
🪸 Annotate Orthogroups
```
Acer_orthogroups_annot <- merge(Acer_orthogroups,Acer_annot_transcripts, by="Transcript")
write.csv(Acer_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Acerhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Mcav_orthogroups_annot <- merge(Mcav_orthogroups,Mcav_annot_transcripts, by="Transcript")
write.csv(Mcav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Mcavhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Ofav_orthogroups_annot <- merge(Ofav_orthogroups,Ofav_annot_transcripts, by="Transcript")
write.csv(Ofav_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Ofavhost_orthogroup_tx2gene_annot.csv",row.names = FALSE)

Past_orthogroups_annot <- merge(Past_orthogroups,Past_annot_transcripts, by="Transcript")
write.csv(Past_orthogroups_annot,file="~/Desktop/Microbial Metatranscriptomics/OrthoFinder/Host/SingleBestORF/Pasthost_orthogroup_tx2gene_annot.csv",row.names = FALSE)


```

</details>
<details>
<summary>BBMAP - Index References for BBsplit</summary>

# BBMAP

> References need to be indexed before bbsplit can be ran. 
```
# Acer
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_Acer_algae_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa

# Mcav 
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_coral_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_Mcav_algae_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_Bacteria_only_transcriptome.fa

# Ofav 
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_coral_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_Ofav_algae_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_Bacteria_only_transcriptome.fa

# Past 
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_coral_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_Past_algae_only_transcriptome.fa
/home/cns.local/nicholas.macknight/software/bbmap/bbmap.sh ref=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_Bacteria_only_transcriptome.fa
```
</details>

<details>
<summary>BBSplit - Assign Reads to Holobiont References</summary>

# BBSplit
### Jan 20th 2025 BBSPLIT

#### BBSplit

1/8/24: After digging into bbsplit methods (which involved a lot of optimization) here is what was performed. 

Reliable Algal symbiont references were concatendated 
Algal symbiont databases were made based on concatenated clade specific references.
Trinity metatranscriptomes were blastn against symbiont databases.
mapped transcripts were quality filtered.
the sequences of QC passed transcripts were extracted via cdbyank from the metatranscritpomes to make clade only references. 
These clade only references were input for bbsplit. 

### End 1/8/24 Notes to self.

### Acropora cervicornis

> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Acer
> Acer_bbsplit_CladeARef.sh
```
#!/bin/bash

# Add directories to the PATH
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/

# Define the reference directory
DIR=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references
HDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer
SDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Acer

# Number of threads to use
threads=50

# Process each FASTQ file
find ${SDIR} -name "*_R1_clean_merged.fastq.gz" | while read FILE; do
    echo "Processing: ${FILE}"
    SAMP=$(basename ${FILE} _R1_clean_merged.fastq.gz)
    echo "Sample: ${SAMP}"

    # Run bbsplit for each sample with increased threads
    /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh \
    in1="${SDIR}/${SAMP}_R1_clean_merged.fastq.gz" \
    in2="${SDIR}/${SAMP}_R2_clean_merged.fastq.gz" \
    ref="${HDIR}/Acer_coral_only_transcriptome.fa,${HDIR}/Acer_Bacteria_only_transcriptome.fa,${DIR}/Clade_A_Acer_algae_only_transcriptome.fa"\
    basename="${SAMP}_%.fq.gz" \
    refstats="${SAMP}_stats.txt" \
    ambig=best ambig2=best \
    outu1="${SAMP}_bboutu_1.fq.gz" \
    outu2="${SAMP}_bboutu_2.fq.gz" \
    threads=${threads}

    # Wait for the current process to finish before moving to the next sample
    wait
done
```

Split Mapped reads from BBsplit into Forward and Reverse reads
> 1_2.sh
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in *.fq.gz; do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=Acer_output_FR/${SAMP}_1.fq out2=Acer_output_FR/${SAMP}_2.fq
done
```

### Montastraea cavernosa

> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Mcav/
> Mcav_bbsplit_CladeCRef.sh
```
#!/bin/bash

# Add directories to the PATH
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/

# Define the reference directory
DIR=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references
HDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav
SDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Mcav

# Number of threads to use
threads=50

# Process each FASTQ file
find ${SDIR} -name "*_R1_clean_merged.fastq.gz" | while read FILE; do
    echo "Processing: ${FILE}"
    SAMP=$(basename ${FILE} _R1_clean_merged.fastq.gz)
    echo "Sample: ${SAMP}"

    # Run bbsplit for each sample with increased threads
    /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh \
    in1="${SDIR}/${SAMP}_R1_clean_merged.fastq.gz" \
    in2="${SDIR}/${SAMP}_R2_clean_merged.fastq.gz" \
    ref="${HDIR}/Mcav_coral_only_transcriptome.fa,${HDIR}/Mcav_Bacteria_only_transcriptome.fa,${DIR}/Clade_C_Mcav_algae_only_transcriptome.fa"\
    basename="${SAMP}_%.fq.gz" \
    refstats="${SAMP}_stats.txt" \
    ambig=best ambig2=best \
    outu1="${SAMP}_bboutu_1.fq.gz" \
    outu2="${SAMP}_bboutu_2.fq.gz" \
    threads=${threads}

    # Wait for the current process to finish before moving to the next sample
    wait
done

```

> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Mcav/
```
mkdir Mcav_output_FR
mv *.fq.gz Mcav_output_FR/
mv stats.txt Mcav_output_FR/
```

Split Mapped reads from BBsplit into Forward and Reverse reads
> 1_2.sh
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in *.fq.gz; do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=Mcav_output_FR/${SAMP}_1.fq out2=Mcav_output_FR/${SAMP}_2.fq
done
```

### Orbicella faveolata

> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Ofav/
> Ofav_bbsplit_CladeDOfavRef.sh
```
#!/bin/bash

# Add directories to the PATH
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/

# Define the reference directory
DIR=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references
HDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output
SDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Ofav

# Number of threads to use
threads=50

# Process each FASTQ file
find ${SDIR} -name "*_R1_clean_merged.fastq.gz" | while read FILE; do
    echo "Processing: ${FILE}"
    SAMP=$(basename ${FILE} _R1_clean_merged.fastq.gz)
    echo "Sample: ${SAMP}"

    # Run bbsplit for each sample with increased threads
    /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh \
    in1="${SDIR}/${SAMP}_R1_clean_merged.fastq.gz" \
    in2="${SDIR}/${SAMP}_R2_clean_merged.fastq.gz" \
    ref="${HDIR}/Ofav_coral_only_transcriptome.fa,${HDIR}/Ofav_Bacteria_only_transcriptome.fa,${DIR}/Clade_D_Ofav_algae_only_transcriptome.fa" \
    basename="${SAMP}_%.fq.gz" \
    refstats="${SAMP}_stats.txt" \
    ambig=best ambig2=best \
    outu1="${SAMP}_bboutu_1.fq.gz" \
    outu2="${SAMP}_bboutu_2.fq.gz" \
    threads=${threads}

    # Wait for the current process to finish before moving to the next sample
    wait
done

```
> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Ofav/
```
mkdir Ofav_output_FR
mv *.fq.gz Ofav_output_FR/
mv stats.txt Ofav_output_FR/
```

Split Mapped reads from BBsplit into Forward and Reverse reads
> 1_2.sh
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in *.fq.gz; do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=Ofav_output_FR/${SAMP}_1.fq out2=Ofav_output_FR/${SAMP}_2.fq
done
```

### Porites astreoides

> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Past/
> Past_bbsplit_CladeCPastRef.sh
```
#!/bin/bash

# Add directories to the PATH
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/

# Define the reference directory
DIR=/home/cns.local/nicholas.macknight/references/Algal_Symbiont_references
HDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past
SDIR=/home/cns.local/nicholas.macknight/SCTLDRNA/MergedFastpProcessedData/Past

# Number of threads to use
threads=50

# Process each FASTQ file
find ${SDIR} -name "*_R1_clean_merged.fastq.gz" | while read FILE; do
    echo "Processing: ${FILE}"
    SAMP=$(basename ${FILE} _R1_clean_merged.fastq.gz)
    echo "Sample: ${SAMP}"

    # Run bbsplit for each sample with increased threads
    /home/cns.local/nicholas.macknight/software/bbmap/bbsplit.sh \
    in1="${SDIR}/${SAMP}_R1_clean_merged.fastq.gz" \
    in2="${SDIR}/${SAMP}_R2_clean_merged.fastq.gz" \
    ref="${HDIR}/Past_coral_only_transcriptome.fa,${HDIR}/Past_Bacteria_only_transcriptome.fa,${DIR}/Clade_A_Past_algae_only_transcriptome.fa"\
    basename="${SAMP}_%.fq.gz" \
    refstats="${SAMP}_stats.txt" \
    ambig=best ambig2=best \
    outu1="${SAMP}_bboutu_1.fq.gz" \
    outu2="${SAMP}_bboutu_2.fq.gz" \
    threads=${threads}

    # Wait for the current process to finish before moving to the next sample
    wait
done
```
> /home/cns.local/nicholas.macknight/SCTLDRNA/bbsplit/Past/
```
mkdir Past_output_FR
mv *.fq.gz Past_output_FR/
mv stats.txt Past_output_FR/
```

Split Mapped reads from BBsplit into Forward and Reverse reads
> 1_2.sh
```
#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/bbmap/
PATH=$PATH:/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/
for FILE in *.fq.gz; do
        echo ${FILE}
        SAMP=$(basename -s .fq.gz $FILE)
        echo $SAMP
/home/cns.local/nicholas.macknight/.sdkman/candidates/java/current/bin/java -ea -Xmx10g  -cp  /home/cns.local/nicholas.macknight/software/bbmap/current/ jgi.ReformatReads in=${SAMP}.fq out1=Past_output_FR/${SAMP}_1.fq out2=Past_output_FR/${SAMP}_2.fq
done
```
</details>

<details>
<summary>Salmon Indexing</summary>

# Salmon

> Index and then perform Quantification

> First have to build a salmon index for your transcriptome. Assume that transcripts.fa contains the set of transcripts you wish to quantify. 

> The purpose of indexing the transcriptomes before running Salmon is to prepare the reference transcriptome in a way that enables fast and efficient quantification of RNA-seq reads. Indexing is a crucial preprocessing step that allows Salmon to rapidly match sequencing reads to the transcriptome without performing traditional alignment, significantly improving speed and computational efficiency.

🪸 Acer Host
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Acer/Acer_index
```

🪸 Mcav Host
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Mcav/Mcav_index
```

🪸 Ofav Host
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav/Ofav_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Ofav/Ofav_index
```

🪸 Past Host
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_coral_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Past/Past_index
```



### Index Building for Symbiont Transcriptomes 
:seedling: Acer Symbiodinium
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_Acer_reference_proteome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Acer/Acer_Clade_A_index
```

:seedling: Mcav Cladicopium 
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_C_Mcav_reference_proteome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Mcav/Mcav_Clade_C_index
```

:seedling: Ofav Durusdinium
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_D_Ofav_algae_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Ofav/Ofav_Clade_D_index
```

:seedling: Past Symbiodinium
```
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/references/Algal_Symbiont_references/Clade_A_Past_algae_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Past/Past_Clade_A_index
```


🦠 Index Building for Bacteria Transcriptomes
```
# Acer
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Acer/Acer_Bacteria_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Bacteria_Acer_index -k 23

# Mcav
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Mcav/Mcav_Bacteria_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Bacteria_Mcav_index -k 23

# Ofav
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Ofav_trinity_output/Ofav_Bacteria_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Bacteria_Ofav_index -k 23

# Past
/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon index -t /home/cns.local/nicholas.macknight/SCTLDRNA/trinity_output_tests/Past/Past_Bacteria_only_transcriptome.fa -i /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/Bacteria_Past_index -k 23
```

count number of transcripts with a count equal or greater to 1
```
awk 'NR>1 && $5 >= 1 { count++ } END { print "Number of transcripts with NumReads >= 1: ", count }' quant.sf
```
</details>
<details>

<summary>Salmon - Quantification</summary>

# Salmon - Quantification

> Salmon processes RNA reads by first matching them to a reference transcriptome. Once the reads are mapped, Salmon counts how many reads align to each gene, providing gene expression counts.

### Acropora cervicornis Salmon Loop
>/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/scripts
>Acer_salmon.sh
```
# Loop for Host

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Acer_coral_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Acer_coral_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Acer_index -l A \
                                -1 ${SAMP}_Acer_coral_only_transcriptome_1.fq \
                                -2 ${SAMP}_Acer_coral_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/${SAMP}_Host_quant
done

# Loop for Algae

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Clade_A_Acer_algae_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Clade_A_Acer_algae_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Acer/Acer_Clade_A_index -l A \
                                -1 ${SAMP}_Clade_A_Acer_algae_only_transcriptome_1.fq \
                                -2 ${SAMP}_Clade_A_Acer_algae_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/${SAMP}_Algae_quant
done
# Loop for Bacteria

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Acer_Bacteria_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Acer_Bacteria_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Bacteria_Acer_index -l A \
                                -1 ${SAMP}_Acer_Bacteria_only_transcriptome_1.fq \
                                -2 ${SAMP}_Acer_Bacteria_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/${SAMP}_Bacteria_quant
done
```

### Montastraea cavernosa Salmon Loop
>/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/scripts
>Mcav_salmon.sh
```
# Loop for Host

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Mcav_coral_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Mcav_coral_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Mcav_index -l A \
				-1 ${SAMP}_Mcav_coral_only_transcriptome_1.fq \
				-2 ${SAMP}_Mcav_coral_only_transcriptome_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/${SAMP}_Host_quant
done

# Loop for Algae

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Clade_C_Mcav_algae_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Clade_C_Mcav_algae_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Mcav/Mcav_Clade_C_index -l A \
				-1 ${SAMP}_Clade_C_Mcav_algae_only_transcriptome_1.fq \
				-2 ${SAMP}_Clade_C_Mcav_algae_only_transcriptome_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/${SAMP}_Algae_quant
done
# Loop for Bacteria

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Mcav_Bacteria_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Mcav_Bacteria_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Bacteria_Mcav_index -l A \
				-1 ${SAMP}_Mcav_Bacteria_only_transcriptome_1.fq \
				-2 ${SAMP}_Mcav_Bacteria_only_transcriptome_2.fq \
				-p 8 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/${SAMP}_Bacteria_quant
done
```

### Orbicella faveolata Salmon Loop
>/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/scripts
>Ofav_salmon.sh
```
# Loop for Host

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Ofav_coral_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Ofav_coral_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Ofav_index -l A \
                                -1 ${SAMP}_Ofav_coral_only_transcriptome_1.fq \
                                -2 ${SAMP}_Ofav_coral_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/${SAMP}_Host_quant
done

# Loop for Algae

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Clade_D_Ofav_algae_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Clade_D_Ofav_algae_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Ofav/Ofav_Clade_D_index -l A \
                                -1 ${SAMP}_Clade_D_Ofav_algae_only_transcriptome_1.fq \
                                -2 ${SAMP}_Clade_D_Ofav_algae_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/${SAMP}_Algae_quant
done
# Loop for Bacteria

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Ofav_Bacteria_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Ofav_Bacteria_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Bacteria_Ofav_index -l A \
                                -1 ${SAMP}_Ofav_Bacteria_only_transcriptome_1.fq \
                                -2 ${SAMP}_Ofav_Bacteria_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/${SAMP}_Bacteria_quant
done
```

### Porites astreoides Salmon Loop
>/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/scripts
>Past_salmon.sh
```
# Loop for Host

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Past_coral_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Past_coral_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Past_index -l A \
                                -1 ${SAMP}_Past_coral_only_transcriptome_1.fq \
                                -2 ${SAMP}_Past_coral_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/${SAMP}_Host_quant
done

# Loop for Algae

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Clade_A_Past_algae_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Clade_A_Past_algae_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Past/Past_Clade_A_index -l A \
                                -1 ${SAMP}_Clade_A_Past_algae_only_transcriptome_1.fq \
                                -2 ${SAMP}_Clade_A_Past_algae_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/${SAMP}_Algae_quant
done
# Loop for Bacteria

#!/bin/bash
PATH=$PATH:/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin
for FILE in ./$(ls *_Past_Bacteria_only_transcriptome_1.fq); do
        echo ${FILE}
        SAMP=$(basename -s _Past_Bacteria_only_transcriptome_1.fq $FILE)
        echo $SAMP
        DIR=/home/cns.local/nicholas.macknight/SCTLDRNA/salmon/

/home/cns.local/nicholas.macknight/software/salmon-1.5.2_linux_x86_64/bin/salmon quant -i ${DIR}/Bacteria_Past_index -l A \
                                -1 ${SAMP}_Past_Bacteria_only_transcriptome_1.fq \
                                -2 ${SAMP}_Past_Bacteria_only_transcriptome_2.fq \
                                -p 50 --validateMappings -o /home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/${SAMP}_Bacteria_quant
done
```









</details>
<details>

<summary>Export</summary>

The primary datasets of interest for exporting are the **salmon quantification results**, **OrthoFinder Results**, and **Annotations**. 
While all the previous code has been conducted within the server, these commands are happening on my local computer to extract these datasets from the server to my lcoal computer.

## Salmon

### Host
```
# Acer
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Acer

# Mcav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Mcav

# Ofav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Ofav

# Past
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/*_Host_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/host/Past

```


### Algal Symbiont
```
# Clade_A_Acer
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/*Algae_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/algae/Clade_A_Acer

# Clade_C_Mcav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/*Algae_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/algae/Clade_C_Mcav


# Clade_D_Ofav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/*Algae_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/algae/Clade_D_Ofav


# Clade_A_Past
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/*Algae_quant /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/algae/Clade_A_Past

```
### Bacteria
```
# Acer
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Acer/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Acer

# Mcav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Mcav/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Mcav

# Ofav
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Ofav/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Ofav

# Past
scp -r nicholas.macknight@holocron:../../home/cns.local/nicholas.macknight/SCTLDRNA/salmon/quants/Past/*Bacteria* /Users/nicholas.macknight/Desktop/Microbial\ Metatranscriptomics/salmon/bacteria/Bacteria_Past
```



</details>
