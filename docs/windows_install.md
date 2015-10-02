####Windows Install Instructions####
___

Note: This has been tested with Win 7 (as that was my only test machine), there maybe other steps required for different windows configurations.

#####1) Need to install Python 2.7#####

Download appropriate installer here:
https://www.python.org/downloads/release/python-2710/


#####2) Install Git for Windows#####
Download here: https://git-scm.com/download/win
* Make sure to make git tools available to windows command line during install

#####3) Open up command line in administrator mode
* Search for 'cmd' from start menu, right click and Run as Administrator

Install Python dependencies:
```
pip install biopython natsort
```

#####4) Now download this repository using git
```
#move into folder of choice
cd C:\Program Files

#clone repository
git clone https://github.com/nextgenusfs/ufits
```

#####5) Download/Install USEARCH8 - get it here: http://www.drive5.com/usearch/download.html.
* copy usearch8 exe file into ufits folder (C:\Program Files\ufits)
* change file name to usearch8 (right click and rename)

#####6) Add location of scripts to PATH variable
Short instructions:
* Open System Control Panel (Start - Settings - Control Panel - System)
* Select Advanced tab
* Open Environmental Variables
* Edit System Variables, Path
* Add to end of string, ;C:\Program Files\ufits
* Now close window
You will need to restart the command prompt for the new settings to work.

#####6) Test Installation
Open command prompt, navigate to the `test_data` folder of ufits.

```
#test scripts on Ion PGM data
ufits-process_ion.py -i ion.test.fastq -o ion
#run clustering
ufits-OTU_cluster.py -i ion.demux.fq -o ion --uchime_ref ITS2 --mock BC_5
```
```
#test scripts on MiSeq data
ufits-process_illumina_folder.py -i illumina_test_data/
#run clustering
ufits-OTU_cluster.py -i ufits.demux.fq -o miseq --uchime_ref ITS2 --mock spike
```



