####Mac OSX Install Instructions####
___

Note: there are many ways to install these dependencies, this is one way.

#####1) Need to install Xcode (Developer Tools that Apple doesn’t install natively but supports):#####
```
xcode-select --install
```

#####2) Install HomeBrew (copy and paste this into terminal):#####
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)”
```

Then setup homebrew: type `brew doctor`, then type: `brew tap homebrew/science`

#####3) Install some tools

```
#install python via homebrew, if not already installed
brew install python

#use pip to install biopython, etc might require sudo
pip install biopython natsort pandas numpy matplotlib
```

#####4) Now install AMPtk via homebrew:
```
#tap my homebrew repository
brew tap nextgenusfs/tap

#install amptk
brew install amptk
```

#####5) You will also need to install USEARCH9 - get it [here](http://www.drive5.com/usearch/download.html).  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch9.2.64_i86osx32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch9.2.64_i86osx32 /usr/local/bin/usearch9
```

#####6) Test Installation
Open terminal, navigate to the `test_data` folder of amptk.  If you installed on Mac with homebrew, it should be here: `/usr/local/opt/amptk/libexec/test_data`.

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


