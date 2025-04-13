# ls split_files/*.fasta | parallel -j 12 "makeblastdb -in {} -dbtype nucl -out ./db/{.}_db -max_file_sz 4GB"
# ls 40_split_files/*.fasta | parallel -j 40 "makeblastdb -in {} -dbtype nucl -out ./db/{.}_db -max_file_sz 4GB -parse_seqids"
seqkit fx2tab -l -n Ldavi.chr.fasta > seq_lengths_per_sequence.txt
python3 generate_bed.py 
seqkit subseq --bed split_sequences.bed Ldavi.chr.fasta > split_Ldavi.chr.fasta
makeblastdb -in split_Ldavi.chr.fasta -dbtype nucl -out lily_db
query="PAIR3-g5032.fa"
# 功能1:
blastn -query $query -db lily_db -out $query.raw -evalue 1e-5 -outfmt 6 -num_threads 10
python3 conversion_blast_results.py $query.raw $query.results
rm $query.raw
# 功能2
seqkit subseq --bed {} -o {}.fa
# 蛋白查询
# makeblastdb -in lily.pep.fa -dbtype prot -out lily_prot_db # 已完成
blastp -query $query -db lily_prot_db -out $query.prot.blast.results -evalue 1e-5 -outfmt 6 -num_threads 10
