###Methods for creating Pre-packaged Databases in UFITS

####Fungal ITS DB
These databases were created from Unite v7.1, first downloading two databases from the UNITE website.  First the General FASTA release of the DB (here)[https://unite.ut.ee/sh_files/sh_general_release_s_22.08.2016.zip] and secondly the Full UNITE+INSD database (here)[https://unite.ut.ee/sh_files/UNITE_public_22.08.2016.fasta.zip].  The taxonomy information is then reformated and databases produced as follows:
```
#Create ITS1 USEARCH Database, trim primers if found, convert taxonomy, and create USEARCH database
ufits database -i UNITE_public_22.08.2016.fasta -o ITS -f ITS1-F -r ITS4 --create_db usearch --keep_all

#Create UTAX Databases
ufits database -i sh_general_release_dynamic_s_22.08.2016_dev.fasta -o ITS_UTAX --create_db utax -f ITS1-F -r ITS4 --keep_all
ufits database -i sh_general_release_dynamic_s_22.08.2016_dev.fasta -o ITS1_UTAX --create_db utax -f ITS1-F -r ITS2 --keep_all
ufits database -i sh_general_release_dynamic_s_22.08.2016_dev.fasta -o ITS2_UTAX --create_db utax -f fITS7 -r ITS4
```

####Arthropod/Chordate mtCOI DB
These data were pulled from the (BOLDv4 database)[http://v4.boldsystems.org].  Since most studies using mtCOI regions are interested in identification of insects in the diets of animals, the BOLD database was queried as follows.  All Chordata sequences were downloaded by querying the (BIN database using the search term Chordata)[http://v4.boldsystems.org/index.php/Public_BINSearch?query=Chordata&searchBIN=Search+BINs].  Similarly, the Arthropods were searched by querying the (BIN databases using the search term Arthropoda)[http://v4.boldsystems.org/index.php/Public_BINSearch?query=Arthropoda&searchBIN=Search+BINs].  All data was then downloaded as TSV output.

The TSV output files (~ 6GB) where then each formatted using the following method, which reformats the taxonomy information and pulls sequences that are annotated in BINS.
```
#reformat taxonomy
ufits/util/bold2utax.py -i Arthropoda_bold_data.txt -o arthropoda.bold.bins.fa
ufits/util/bold2utax.py -i Chordata_bold_data.txt -o chordata.bold.bins.fa

#combine datasets
cat arthropoda.bold.bins.fa chordata.bold.bins.fa > all.data.bins.fa
```
The data is then processed with a second script that will loop through each valid BIN and pull out the centroid cluster of the BIN using UCLUST alignment.  The script will then trim the forward primer if it exists and truncate the data to the first 200 bp of the COI region - this is the region that is typically targeted with published datasets.  The script will then output two files, 1 for use in creation of a USEARCH database and the other for creation of a UTAX database.
```
ufits/util/bold2ufits.py -i all.data.bins.fa -o arthropods.chordates
```
Finally, the data is processed for UFITS database command:
```
ufits database -i arthropods.chordates.all4usearch.fa -o COI --create_db usearch --keep_all --skip_trimming --format off
ufits database -i arthropods.chordates.all4utax.fa -o COI_UTAX --create_db utax --keep_all --skip_trimming --format off
```

####16S database
This is downloaded from R. Edgar's (website)[http://drive5.com/utax/data/rdp_v16.tar.gz] and then formatted for UFITS:
```
ufits database -i rdp_v16.fa -o 16S --format off --create_db utax --skip_trimming --keep_all
```


