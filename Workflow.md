# Workflow

The objectives of this study are to assess the distribution, diversity, abundance, and phylogenetic distribution of 1-Aminocyclopropane-1-Carboxylate Deaminase (ACC-deamianse) in the envrionment. 

This workflow was performed with QIIME, MUSCLE, and raxML on a Liunux-based machine.

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

