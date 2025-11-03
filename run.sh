# ls split_files/*.fasta | parallel -j 12 "makeblastdb -in {} -dbtype nucl -out ./db/{.}_db -max_file_sz 4GB"
# ls 40_split_files/*.fasta | parallel -j 40 "makeblastdb -in {} -dbtype nucl -out ./db/{.}_db -max_file_sz 4GB -parse_seqids"
seqkit fx2tab -l -n Ldavi.chr.fasta > seq_lengths_per_sequence.txt
python3 generate_bed.py 
seqkit subseq --bed split_sequences.bed Ldavi.chr.fasta > split_Ldavi.chr.fasta
makeblastdb -in split_Ldavi.chr.fasta -dbtype nucl -out lily_db
makeblastdb -in Ldavi.chr.pep.fa -dbtype prot -out lily_prot_db

query="REC8cds.fasta"
# 功能1:
blastn -query $query -db lily_db -out $query.raw.txt -evalue 1e-5 -outfmt 6 -num_threads 10
python3 conversion_blast_results.py $query.raw.txt $query.results.txt
rm $query.raw.txt
# 功能2： 序列提取
seqkit subseq --bed uesr.bed Ldavi.chr.fasta > user.fa
# 功能3: 蛋白查询
blastp -query $query -db lily_prot_db -out $query.prot.blast.raw.txt -evalue 1e-5 -outfmt 6 -num_threads 10
python3 conversion_blast_results.py $query.prot.blast.raw.txt $query.prot.blast.results.txt