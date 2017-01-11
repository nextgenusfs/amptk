####Linux Install Instructions####
___

Note: there are many ways to install these dependencies, this is one way.

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

#####5) Test Installation
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


