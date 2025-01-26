pwd
# /config/workspace

mkdir -p class-variantcalling

cd class-variantcalling

mkdir -p analysis

cd analysis

mkdir -p raw_data

cd /config/workspace/datiesame

tar -xzvf data_resequencing.tar.gz -C /config/workspace/class-variantcalling/analysis/raw_data

#per ottenere il path corretto di raw_data cliccare con tasto destro su raw_data
#al posto del comando sopra riportato possiamo fare solo tar -xzvf data_resequencing.tar.gz e poi spostare manualmente i campioni

cd ..

cd class-variantcalling/analysis

mkdir -p alignment

cd alignment

## now we can perform the alignment with BWA
##i nomi dei campioni, dopo raw_data, potrebbero cambiare all'esame

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/class-variantcalling/analysis/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/config/workspace/class-variantcalling/analysis/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 2 -bhS -o normal.bam -

## Real time: 176.099 sec; CPU: 256.669 sec

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/class-variantcalling/analysis/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/config/workspace/class-variantcalling/analysis/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 2 -bhS -o disease.bam -

## Real time: 173.232 sec; CPU: 256.204 sec


# sort the bam file
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam

# index the bam file, cotrollare sempre che si siano formate le cartelle con l'estensione.bai
samtools index normal_sorted.bam
samtools index disease_sorted.bam


# Marking duplicates

gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam


### recalibrating

gatk BaseRecalibrator \
   -I normal_md.bam \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/datiesame/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/datiesame/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O normal_recal_data.table

gatk BaseRecalibrator \
   -I disease_md.bam \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/datiesame/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/datiesame/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O disease_recal_data.table


#### Apply recalibration

gatk ApplyBQSR \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I normal_md.bam \
   --bqsr-recal-file normal_recal_data.table \
   -O normal_recal.bam

gatk ApplyBQSR \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I disease_md.bam \
   --bqsr-recal-file disease_recal_data.table \
   -O disease_recal.bam


### variant calling

cd ..

mkdir -p variants

cd variants


## first single sample discovery

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/class-variantcalling/analysis/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/class-variantcalling/analysis/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF

## then consolidate the 2 files

mkdir -p tmp

### on AMD64 this code ######
## combine the files into one
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V normal.g.vcf.gz \
      -V disease.g.vcf.gz \
      --genomicsdb-workspace-path compared_db \
      --tmp-dir /config/workspace/class-variantcalling/analysis/variants/tmp \
      -L chr21

### on AMD64 this code ######
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/datiesame/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /config/workspace/datiesame/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz


#### ANNOTATE THE SAMPLE

#fare pwd per essere sicuri a questo punto di essere nella cartella variants


### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -dataDir /config/workspace/snpeff_data -v hg38 results.vcf.gz >results_ann.vcf


### filter variants

cat results_ann.vcf | grep "#CHROM" | cut -f 10-

grep "#" results_ann.vcf >filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/1\/1/){print $_;}' >>filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/0\/1/){print $_;}' >>filtered_variants.vcf

sudo conda install bioconda::snpsift

#la password per bioconda è student

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "GEN[disease].GT" "GEN[normal].GT" ANN[*].GENE ANN[*].EFFECT

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "REF" "ALT" "GEN[*].GT" ANN[0].GENE ANN[0].EFFECT
