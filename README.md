## Prepare
First download data from cngb. You can find info in https://db.cngb.org/data_resources/project/CNP0005511
```sh
git clone https://github.com/GooLey1025/Lily_Web.git
cd Lily_Web
wget -c ftp://ftp2.cngb.org/pub/CNSA/data5/CNP0005511/CNS1065710/CNA0139751/Ldavi.chr.fasta.gz
wget -c ftp://ftp2.cngb.org/pub/CNSA/data5/CNP0005511/CNS1065710/CNA0139751/Ldavi.chr.pep.fa
bgzip -d -c Ldavi.chr.fasta.gz Ldavi.chr.fasta.gz > Ldavi.chr.fasta
```
Also make sure `seqkit` and `blast` have been installed and added in your environmental variables.

Due to large chromosome length which cause error with `Blast`, we need to split large chromosome into 1GB sub-chromosome.
```sh
seqkit fx2tab -l -n Ldavi.chr.fasta > seq_lengths_per_sequence.txt
python3 generate_bed.py
seqkit subseq --bed split_sequences.bed Ldavi.chr.fasta > split_Ldavi.chr.fasta
# prepare database for blast
makeblastdb -in split_Ldavi.chr.fasta -dbtype nucl -out lily_db
makeblastdb -in Ldavi.chr.pep.fa -dbtype prot -out lily_prot_db
```

## Depoly the web
```sh
pip install -r requirements.txt
# start web server
uvicorn app.main:app --host 0.0.0.0 --port 1111 --reload
```
