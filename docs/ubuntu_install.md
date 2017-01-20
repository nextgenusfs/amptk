####Linux Install Instructions####
___

Note: there are many ways to install these dependencies, I will give you two ways.

####Install using LinuxBrew (might be easier?)
#####1) Install LinuxBrew (copy and paste this into terminal):#####
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
```
Then setup linuxBrew by adding variables to ~/.bashrc or ~/.bash_aliases
```
export PATH="$HOME/.linuxbrew/bin:$PATH"
export MANPATH="$HOME/.linuxbrew/share/man:$MANPATH"
export INFOPATH="$HOME/.linuxbrew/share/info:$INFOPATH"
```
Final setup:  type `brew doctor`.

#####2) Install some necessary tools (may not need to do this depending on your system)
```
#install python via homebrew, if not already installed
brew install python

#use pip to install biopython, etc might require sudo.
pip install biopython natsort pandas numpy matplotlib
```

#####3) Now install AMPtk via linuxbrew:
```
#tap homebrew science
brew tap homebrew/science

#tap my homebrew repository
brew tap nextgenusfs/tap

#install amptk
brew install amptk
```

#####4) You will also need to install USEARCH9 - get it [here](http://www.drive5.com/usearch/download.html).  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch9.2.64_i86osx32
```

```
#create softlink to folder in $PATH, i.e.
sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9
```


####Manual Install on Ubuntu
#####1) Install Python/modules:#####
You should be able to use system python, however, sometimes `pip` can be problematic and can occastionally require `sudo` which some users don't have.
```
#use pip to install biopython, etc might require sudo.
pip install biopython natsort pandas numpy matplotlib biom-format psutil
```
Alternativly you could use Minconda python, which should make installing packages easier.
```
#get the minconda package
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh

#install
bash Miniconda2-latest-Linux-x86_64.sh

#restart terminal window for changes to take place
```
Now you can install python packages with conda
```
conda install biopython natsort pandas numpy matplotlib biom-format psutil
```
#####2) Install External Dependencies:#####
You will need to then get VSEARCH, USEARCH, BEDTOOLS, R, DADA2
```
#install bedtools
sudo apt-get install bedtools

#install vsearch
wget https://github.com/torognes/vsearch/releases/download/v2.3.4/vsearch-2.3.4-linux-x86_64.tar.gz
tar -xvzf vsearch-2.3.4-linux-x86_64.tar.gz
sudo ln -s /path/to/vsearch-2.3.4-linux-x86_64/bin/vsearch /usr/local/bin/vsearch

#install R and DADA2
sudo apt-get install r-base

#start an R session and install DADA2
R

#in the R session do the following
source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates = FALSE)
biocLite("ShortRead", suppressUpdates = FALSE)
biocLite("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2")

#exit R
q()

#Alternativly the amptk dada2 script will try to install DADA2 automatically, however it will not be multithreaded...
```
You will also need to install USEARCH9 - get it [here](http://www.drive5.com/usearch/download.html).  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch9.2.64_i86osx32

#create softlink to folder in $PATH, i.e.
sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9
```
#####3) Install AMPtk:#####
You can now install the amplicon tool kit using git:
```
#clone the repositry
git clone https://github.com/nextgenusfs/amptk.git

#softlink into PATH
sudo ln -s /path/to/amptk/amptk /usr/local/bin/amptk
```

####Test Installation
Open terminal, navigate to the `test_data` folder of amptk. If you installed with LinuxBrew it should be here: `$HOME/.linuxbrew/opt/amptk/libexec/test_data`.

```
#test scripts on Ion PGM data
amptk ion -i ion.test.fastq -o ion

#run clustering
amptk cluster -i ion.demux.fq -o ion
```
```
#test scripts on MiSeq data
amptk illumina -i illumina_test_data/

#run clustering
amptk cluster -i amptk.demux.fq -o miseq
```