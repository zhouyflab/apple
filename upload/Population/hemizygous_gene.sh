###############################
# 1.Build reference indexes
###############################
samtools faidx ./fa/genome.hap1.fa
minimap2 -d genome.hap1.mmi ./fa/genome.hap1.fa
###############################
# 2.Align HiFi reads
###############################
minimap2 -t 30 --cs --MD -ax map-hifi ./fa/genome.hap1.fa ./HiFi/genome.fastq.gz > genome.hap1output.sam
samtools view -bt ./fa/genome.hap1.fa.fai genome.hap1output.sam | samtools sort -@ 30 -o genome.hap1.hifi.sorted.bam -
samtools index -@ 30 genome.hap1.hifi.sorted.bam
###############################
# 3.Detect structural variants
###############################
sniffles     --input genome.hap1.hifi.sorted.bam     --vcf genome.hap1.sv.vcf     --threads 30     --minsupport 4     --minsvlen 50     --reference ./fa/genome.hap1.fa
###############################
# 4.Screening of SV detected by Sniffles
###############################
bcftools view -i 'GT="0/1" && (INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && ABS(INFO/SVLEN)<=6000000 && FILTER="PASS" && INFO/PRECISE=1' genome.hap1.sv.vcf -Ov -o genome.hap1.DEL_DUP_INV.filtered.vcf
###############################
# 5.Identify genes overlapping SV regions
###############################
bcftools view genome.hap1.DEL_DUP_INV.filtered.vcf -Ou | bcftools query -f '%CHROM\t%POS0\t%END\n' | bedtools intersect -a genome.hap1.genes.bed -b stdin -u -f 0.8 > genome.hap1.DEL_DUP_INV.filtered.sv.vcf.hemi.bed
