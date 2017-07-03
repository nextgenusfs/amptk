### Methods for creating Pre-packaged Databases in AMPtk

#### Fungal ITS DB
These databases were created from Unite v7.1.2 (released November 20th, 2016), first downloading two databases from the UNITE website.  First the General FASTA release of the DB [here](https://unite.ut.ee/sh_files/sh_general_release_s_20.11.2016.zip) and secondly the Full UNITE+INSD database [here](https://unite.ut.ee/sh_files/UNITE_public_20.11.2016.fasta.zip).  The taxonomy information is then reformated and databases produced as follows:
```
#Create full length ITS USEARCH Database, trim primers if found, convert taxonomy, and create USEARCH database
amptk database -i UNITE_public_20.11.2016.fasta -o ITS -f ITS1-F -r ITS4 --create_db usearch --keep_all

#Create UTAX Databases
amptk database -i sh_general_release_dynamic_s_20.11.2016_dev.fasta -o ITS_UTAX --create_db utax -f ITS1-F -r ITS4 --keep_all
amptk database -i sh_general_release_dynamic_s_20.11.2016_dev.fasta -o ITS1_UTAX --create_db utax -f ITS1-F -r ITS2 --keep_all
amptk database -i sh_general_release_dynamic_s_20.11.2016_dev.fasta -o ITS2_UTAX --create_db utax -f fITS7 -r ITS4
```

#### Arthropod/Chordate mtCOI DB
These data were pulled from the [BOLDv4 database](http://v4.boldsystems.org)  Since most studies using mtCOI regions are interested in identification of insects in the diets of animals, the BOLD database was queried as follows.  All Chordata sequences were downloaded by querying the [BIN database using the search term Chordata](http://v4.boldsystems.org/index.php/Public_BINSearch?query=Chordata&searchBIN=Search+BINs).  Similarly, the Arthropods were searched by querying the [BIN databases using the search term Arthropoda](http://v4.boldsystems.org/index.php/Public_BINSearch?query=Arthropoda&searchBIN=Search+BINs).  All data was then downloaded as TSV output.

The TSV output files (~ 6GB) where then each formatted using the following method, which reformats the taxonomy information and pulls sequences that are annotated in BINS.
```
#reformat taxonomy
amptk/util/bold2utax.py -i Arthropoda_bold_data.txt -o arthropoda.bold.bins.fa
amptk/util/bold2utax.py -i Chordata_bold_data.txt -o chordata.bold.bins.fa

#combine datasets
cat arthropoda.bold.bins.fa chordata.bold.bins.fa > all.data.bins.fa
```
The data is then processed with a second script that will loop through each valid BIN and pull out the centroid cluster of the BIN using UCLUST alignment.  The script will then trim the forward primer if it exists and truncate the data to the first 200 bp of the COI region - this is the region that is typically targeted with published datasets.  The script will then output two files, 1 for use in creation of a USEARCH database and the other for creation of a UTAX database.
```
amptk/util/bold2amptk.py -i all.data.bins.fa -o arthropods.chordates
```
Finally, the data is processed for AMPtk database command:
```
amptk database -i arthropods.chordates.all4usearch.fa -o COI --create_db usearch --keep_all --skip_trimming --format off
amptk database -i arthropods.chordates.all4utax.fa -o COI_UTAX --create_db utax --keep_all --skip_trimming --format off
```

#### LSU database
The fungal 28S database (LSU) was downloaded from [here](http://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz).  The sequences were then converted into AMPtk databases as follows:
```
amptk database -i fungi.unaligned.fa -o LSU --format rdp2utax --skip_trimming --create_db usearch --derep_fulllength --keep_all
```
To generate a training set for UTAX, the sequences were first dereplicated, and clustered at 97% to get representative sequences for training.  This training set was then converted to a UTAX database:
```
amptk database -i fungi.trimmed.fa -o LSU_UTAX --format off --skip_trimming --create_db utax --keep_all
```


#### 16S database
This is downloaded from R. Edgar's [website](http://drive5.com/utax/data/rdp_v16.tar.gz) and then formatted for AMPtk.  Note there is room for substantial improvement here, I just don't typically work on 16S - so please let me know if you want some suggestions on what to do here.
```
amptk database -i rdp_v16.fa -o 16S --format off --create_db utax --skip_trimming --keep_all
```


