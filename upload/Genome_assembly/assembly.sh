###### primiary_assemble
hifiasm -o hifiasm_out -t10 --h1 $hic1 --h2 $hic2 --ul $ont $hifi

###### HiC_scaffolding
bwa index genome.hap1hap2.fa
python ~/software/juicer-1.6/misc/generate_site_positions.py DpnII genome genome.hap1hap2.fa 
awk 'BEGIN{OFS="\t"}{print $1, $NF}' genome.hap1hap2_DpnII.txt > .genome.hap1hap2.chrom.sizes
~/software/juicer-1.6/CPU/juicer.sh -g hap1_hic -z genome.hap1hap2.fa -s DpnII -y genome.hap1hap2_DpnII.txt -p genome.hap1hap2.chrom.sizes -t 56 
~/software/3d-dna/run-asm-pipeline.sh  -r 0  reference/genome.hap1hap2.fa  aligned/merged_nodups.txt
~/software/3d-dna/run-asm-pipeline-post-review.sh -r genome.hap1hap2.review.assembly ./reference/genome.hap1hap2.fa aligned/merged_nodups.txt



