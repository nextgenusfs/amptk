####Ubuntu Install Instructions####
___

Note: there are many ways to install these dependencies, this is one way.

#####1) Install Dependencies

```
#make sure python is installed
sudo apt-get install python

#use pip to install biopython and natsort, might require sudo
pip install biopython natsort  
```

#####2) Now download this repository:

`git clone https://github.com/nextgenusfs/ufits`

And then you will need to add to your `~/.bash_aliases` or always include the entire path to the scripts at runtime.

`export PATH="/location/of/packages/ufits:$PATH"`

#####3) You will also need to install USEARCH8 - get it here: http://www.drive5.com/usearch/download.html.  One way to make the program executable and move into your path:

```
#make executable
sudo chmod +x /path/to/usearch8.0.1623_i86linux32
```

```
#create softlink to folder in $PATH
sudo ln -s /path/to/usearch8.0.1623_i86linux32 /usr/local/bin/usearch8

#or rename it and move to folder in $PATH
sudo mv /path/to/usearch8.0.1623_i86linux32 /usr/local/bin/usearch8
```

#####4) Test Installation
Open command prompt, navigate to the `test_data` folder of ufits.

```
#test scripts on Ion PGM data
ufits ion -i ion.test.fastq -o ion
#run clustering
ufits cluster -i ion.demux.fq -o ion --uchime_ref ITS2 --mock BC_5
```
```
#test scripts on MiSeq data
ufits illumina -i illumina_test_data/
#run clustering
ufits cluster -i ufits.demux.fq -o miseq --uchime_ref ITS2 --mock spike
```

