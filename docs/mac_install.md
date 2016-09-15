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

#####4) Now install UFITS via homebrew:
```
#tap my homebrew repository
brew tap nextgenusfs/tap

#install ufits
brew install ufits
```

#####5) You will also need to install USEARCH8 - get it [here](http://www.drive5.com/usearch/download.html).  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch8.1.1861_i86osx32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch8.1.1861_i86osx32 /usr/local/bin/usearch8
```

#####6) Test Installation
Open terminal, navigate to the `test_data` folder of ufits.

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


