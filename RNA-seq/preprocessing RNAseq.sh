#checks quality of fastq files
fastqc *.gz -t 4 #number of threads it's using, can be changed
echo "fastqc finished"

#builds index for reference genome
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa #change name if a different version, do the same downstream
echo "Index built"

#pseudoalignment and output in separate directories
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A10 -t 4 A10_1.fq.gz A10_2.fq.gz &> A10.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A11 -t 4 A11_1.fq.gz A11_2.fq.gz &> A11.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A12 -t 4 A12_1.fq.gz A12_2.fq.gz &> A12.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A13 -t 4 A13_1.fq.gz A13_2.fq.gz &> A13.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A14 -t 4 A14_1.fq.gz A14_2.fq.gz &> A14.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A15 -t 4 A15_1.fq.gz A15_2.fq.gz &> A15.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A16 -t 4 A16_1.fq.gz A16_2.fq.gz &> A16.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A17 -t 4 A17_1.fq.gz A17_2.fq.gz &> A17.log
kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o A18 -t 4 A18_1.fq.gz A18_2.fq.gz &> A18.log


echo "Alignment finished"

#summarizes qc and mapping results in single html
multiqc -d .