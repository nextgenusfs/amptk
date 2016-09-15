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

#####3) Now install UFITS via linuxbrew:
```
#tap homebrew science
brew tap homebrew/science

#tap my homebrew repository
brew tap nextgenusfs/tap

#install ufits
brew install ufits
```

#####4) You will also need to install USEARCH8 - get it [here](http://www.drive5.com/usearch/download.html).  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch8.1.1861_i86linux32
```

```
#create softlink to folder in $PATH, i.e.
sudo ln -s /path/to/usearch8.1.1861_i86osx32 /usr/local/bin/usearch8
```

#####5) Test Installation
Open terminal, navigate to the `test_data` folder of ufits. If you installed with LinuxBrew it should be here: `$HOME/.linuxbrew/opt/ufits/libexec/test_data`.

```
#test scripts on Ion PGM data
ufits ion -i ion.test.fastq -o ion

#run clustering
ufits cluster -i ion.demux.fq -o ion
```
```
#test scripts on MiSeq data
ufits illumina -i illumina_test_data/

#run clustering
ufits cluster -i ufits.demux.fq -o miseq
```


