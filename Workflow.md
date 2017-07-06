# Workflow

The objectives of this study are to assess the distribution, diversity, abundance, and phylogenetic distribution of 1-Aminocyclopropane-1-Carboxylate Deaminase (ACC-deamianse) in the envrionment. 

This workflow was performed with QIIME, MUSCLE, R, and raxML on a Liunux-based machine.

## Getting data

Sequences were aggregated from JGI's IMG M server (https://img.jgi.doe.gov/cgi-bin/m/main.cgi) which aggregates data sequenced at JGI as well as from NCBI.

Seqeunces were gathered by using the gene search function, searching for "1-Aminocyclopropane-1-Carboxylate Deaminase" and pulling 50 studies worth of data at a time. Metadata for each sequence pulled was also gathered.

In addition, using the KEGG pathways option, I gathered abundance of ACC-deaminase in metagenomes as well as the single copy gene rplB to normalize ACC-deaminase abundance accross the studies. 

This leaves us with 4 files:

        1) full_nulc.fasta- ~150K ACC sequences
        2) ACC_rplB_comb.txt- file with the abundance of ACC and rplB in metagenomes.
        3) collated_metadata.txt- the metadata for all 150K sequences
        4) collapsed_metadata.txt- the metadata for the studies data was pulled from, will use as a QIIME mapping file

## Scripts

First thing we need to do is purge the dataset of non-ACC deaminase sequences. These are seqeuences like "hypothetical protein" or "ACC cysteine residue". These aren't true ACC deaminase and will affect our analysis. 

We'll use QIIME's 'filter_fasta.py' command to do this paired with a list of bad sequence identifires (can get from metadata file).

```
filter_fasta.py -n -i full_nulc.fasta -s not_acc_seqs.txt -o filted_nulc.fasta
```

If we do grep we can see we lost about ~7k sequences.

```
grep -c ">" full_nulc.fasta
157,730

grep -c ">" filted_nulc.fasta
150,831
```
Next we'll collapse the sequences down into OTUs at a 85% sequence identity using swarm (Mahe et al. 2013, PeerJ).
```
pick_otus.py -s 0.85 -m swarm -i filted_nulc.fasta -o 85_otus
```

Next we'll pull a set of representative sequence for each of the OTUs.
```
pick_rep_set.py -i 85_otus/filted_nulc_otus.txt -o acc_rep_set.fna -f filted_nulc.fasta
grep -c ">" acc_rep_set.fna
26,297
```

Make an OTU table, calculate alpha diversity, and append the alpha diversity to the mapping file for easier analysis in R. 
```
make_otu_table.py -i 85_otus/filted_nulc_otus.txt -o acc_otu_table.biom
alpha_diversity.py -i acc_otu_table.biom -o acc_alpha.txt -m observed_species,shannon,singles
add_alpha_to_mapping_file.py -m acc_map.txt -i acc_alpha.txt -o acc_alphaphied_map.txt
```

Next we want to use BLASTn to assign our OTUs taxonomy using the nucleotide database
```
#first, get the most up-to-date nt database
update_blastdb.pl nt

#unzip the files
tar -xzvf nt*

blastn -query acc_rep_set.fna -max_target_seqs 1 -outfmt "6 qseqid sacc stitle pident evalue" -out acc_results_out -negative_gilist sequence.gi -db nt

#sequence.gi is a file of envrionmental sequences that was retrieved from NCBI Entrez:

https://www.ncbi.nlm.nih.gov/nuccore/?term=%22environmental+samples%22%5Borganism%5D+OR+metagenomes%5Borgn%5D+OR+%22Unidentified%22+OR+%22clone%22
```
This produces only a list of ascession numbers and the name that corresponds to the ascession number. We want to map this data back to a taxonomic assignment. The details of which can be found in my GitHub repo (BLAST-taxonomy). 


## Making figures & Statistics

We'll first start off by making a map that has all the studies plotted on a map of the earth. 

```
library(ggplot2)
library(readr)
library(rworldmap)

#read in data
acc_map <- read_delim("~/acc_map.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
str(acc_map)

#make columns of lat/long numeric
acc_map$lat<-as.numeric(acc_map$lat)
acc_map$long<-as.numeric(acc_map$long)

#get world map data
worldmap <- map_data("world")
str(worldmap)

# Plot entire world map with studies with ACC
  ggplot(data = worldmap, aes(x = long, y = lat, group = group)) +
  geom_polygon(colour = "grey", fill = "grey") +
  theme_bw() +
  geom_point(data = acc_map, aes(x = long, y = lat, fill=Ecosystem_Type, colour = Ecosystem_Type, group=Ecosystem_Type), size = 3)
```

Next we'll make a box and wisker plot of the abundance of ACC deaminase accross different environments.

```
library(readr)
library(ggplot2)

acc_abun <- read_delim("~/acc_abun.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

ggplot(acc_abun, aes(Order, Normalied_ACC, fill=Order))+
  geom_boxplot()+
  ylab("Normailzed ACC Abundance")+
  xlab("")+
  theme_bw()+
  coord_cartesian(ylim=c(0,1))+
  coord_flip()
```

Make a box and wisker plot of the alpha diversity of ACC deaminase in different environments.

```
library(readr)
library(ggplot2)

#read data
acc_alphaphied_map2 <- read_delim("~/acc_alphaphied_map2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#plot OTUs observed
ggplot(acc_alphaphied_map2, aes(Ecosystem_Type, otus_obs, fill=Ecosystem_Type))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  scale_y_log10()

#plot singeltons
ggplot(acc_alphaphied_map2, aes(Ecosystem_Type, singles, fill=Ecosystem_Type))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  scale_y_log10()

#plot Shannon Diversity
ggplot(acc_alphaphied_map2, aes(Ecosystem_Type, shannon, fill=Ecosystem_Type))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()+
  scale_y_log10()
```
