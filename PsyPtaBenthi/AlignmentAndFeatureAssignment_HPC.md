# Pseudomonas-N. benthamiana RNA-Seq Pipeline
***Morgan Carter 2021-2022***

## HISAT2 v2.2.1 (HPC) - alignment for N. benthamiana
Using [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) for alignment. First need to make an index of the various genomes using the `hisat2-build` command. It's just an input and an output. 

For the input, I'm using the [Nb101 genome](https://solgenomics.net/organism/Nicotiana_benthamiana/genome) and using the re-annotated version (Additional File 5) from [this paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6058-6#Tab2) :
Kourelis, J., Kaschani, F., Grosse-Holz, F.M. _et al._ A homology-guided, genome-based proteome for improved proteomics in the alloploid _Nicotiana benthamiana_. _BMC Genomics_  **20,** 722 (2019). https://doi.org/10.1186/s12864-019-6058-6

`wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-019-6058-6/MediaObjects/12864_2019_6058_MOESM5_ESM.gff`
`wget ftp://ftp.solgenomics.net/genomes/Nicotiana_benthamiana/assemblies/Niben.genome.v1.0.1.ds.nrcontigs.fasta.gz`

`make_indices.slurm` run using `sbatch`
```bash
#!/bin/bash
#SBATCH --job-name=HISAT2_build
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=01:00:00
module load hisat2
hisat2-build -p 8 /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101
```

Then I did the alignment using the normal `hisat2` command for each of the N. benthamiana samples using that index as an input. The output is a sam file, which is rather large. Immediately I used [samtools](http://samtools.sourceforge.net/) to convert the giant sam file to a bam file (and sort it) and then get rid of the sam file to save space. I also had to do some of this on xdisk space because the lab storage area didn't have enough space. I set the max intron length to 10000 just to try to cut down on crazy long weird alignments.  
`HISAT2_align_nb_B728a.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=HISAT2_align
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH -o %j_HISAT_nb101_B728a.out
#SBATCH --time=03:00:00
module load hisat2
module load samtools
cd /groups/baltrus/morgan/Ptab_B728a_RNA/Nb_B728a_reads/
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-1-P0_S2_L001_R1_001.fastq.gz,Read1/B728a-1-P0_S2_L002_R1_001.fastq.gz,Read1/B728a-1-P0_S2_L003_R1_001.fastq.gz,Read1/B728a-1-P0_S2_L004_R1_001.fastq.gz -2 Read2/B728a-1-P0_S2_L001_R2_001.fastq.gz,Read2/B728a-1-P0_S2_L002_R2_001.fastq.gz,Read2/B728a-1-P0_S2_L003_R2_001.fastq.gz,Read2/B728a-1-P0_S2_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P0_1.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-1-P5_S3_L001_R1_001.fastq.gz,Read1/B728a-1-P5_S3_L002_R1_001.fastq.gz,Read1/B728a-1-P5_S3_L003_R1_001.fastq.gz,Read1/B728a-1-P5_S3_L004_R1_001.fastq.gz -2 Read2/B728a-1-P5_S3_L001_R2_001.fastq.gz,Read2/B728a-1-P5_S3_L002_R2_001.fastq.gz,Read2/B728a-1-P5_S3_L003_R2_001.fastq.gz,Read2/B728a-1-P5_S3_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P5_1.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-2-P0_S6_L001_R1_001.fastq.gz,Read1/B728a-2-P0_S6_L002_R1_001.fastq.gz,Read1/B728a-2-P0_S6_L003_R1_001.fastq.gz,Read1/B728a-2-P0_S6_L004_R1_001.fastq.gz -2 Read2/B728a-2-P0_S6_L001_R2_001.fastq.gz,Read2/B728a-2-P0_S6_L002_R2_001.fastq.gz,Read2/B728a-2-P0_S6_L003_R2_001.fastq.gz,Read2/B728a-2-P0_S6_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P0_2.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-2-P5_S7_L001_R1_001.fastq.gz,Read1/B728a-2-P5_S7_L002_R1_001.fastq.gz,Read1/B728a-2-P5_S7_L003_R1_001.fastq.gz,Read1/B728a-2-P5_S7_L004_R1_001.fastq.gz -2 Read2/B728a-2-P5_S7_L001_R2_001.fastq.gz,Read2/B728a-2-P5_S7_L002_R2_001.fastq.gz,Read2/B728a-2-P5_S7_L003_R2_001.fastq.gz,Read2/B728a-2-P5_S7_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P5_2.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-3-P0_S10_L001_R1_001.fastq.gz,Read1/B728a-3-P0_S10_L002_R1_001.fastq.gz,Read1/B728a-3-P0_S10_L003_R1_001.fastq.gz,Read1/B728a-3-P0_S10_L004_R1_001.fastq.gz -2 Read2/B728a-3-P0_S10_L001_R2_001.fastq.gz,Read2/B728a-3-P0_S10_L002_R2_001.fastq.gz,Read2/B728a-3-P0_S10_L003_R2_001.fastq.gz,Read2/B728a-3-P0_S10_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P0_3.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/B728a-3-P5_S11_L001_R1_001.fastq.gz,Read1/B728a-3-P5_S11_L002_R1_001.fastq.gz,Read1/B728a-3-P5_S11_L003_R1_001.fastq.gz,Read1/B728a-3-P5_S11_L004_R1_001.fastq.gz -2 Read2/B728a-3-P5_S11_L001_R2_001.fastq.gz,Read2/B728a-3-P5_S11_L002_R2_001.fastq.gz,Read2/B728a-3-P5_S11_L003_R2_001.fastq.gz,Read2/B728a-3-P5_S11_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/B728a_nb101/B728a_P5_3.sam
cd /xdisk/baltrus/morgancarter
for i in HISAT2/B728a_nb101/*.sam; do samtools view -bS -@ 16 $i | samtools sort -@ 16 -o $i.sorted.bam && rm $i; done
```
`HISAT2_align_nb_Ptab.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=HISAT2_align
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH -o %j_HISAT_nb101_Ptab.out
#SBATCH --time=03:00:00
module load hisat2
cd /groups/baltrus/morgan/Ptab_B728a_RNA/Nb_Ptab_reads/
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-5-P0_S1_L001_R1_001.fastq.gz,Read1/Ptab-5-P0_S1_L002_R1_001.fastq.gz,Read1/Ptab-5-P0_S1_L003_R1_001.fastq.gz,Read1/Ptab-5-P0_S1_L004_R1_001.fastq.gz -2 Read2/Ptab-5-P0_S1_L001_R2_001.fastq.gz,Read2/Ptab-5-P0_S1_L002_R2_001.fastq.gz,Read2/Ptab-5-P0_S1_L003_R2_001.fastq.gz,Read2/Ptab-5-P0_S1_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P0_1.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-5-P5_S2_L001_R1_001.fastq.gz,Read1/Ptab-5-P5_S2_L002_R1_001.fastq.gz,Read1/Ptab-5-P5_S2_L003_R1_001.fastq.gz,Read1/Ptab-5-P5_S2_L004_R1_001.fastq.gz -2 Read2/Ptab-5-P5_S2_L001_R2_001.fastq.gz,Read2/Ptab-5-P5_S2_L002_R2_001.fastq.gz,Read2/Ptab-5-P5_S2_L003_R2_001.fastq.gz,Read2/Ptab-5-P5_S2_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P5_1.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-6-P0_S3_L001_R1_001.fastq.gz,Read1/Ptab-6-P0_S3_L002_R1_001.fastq.gz,Read1/Ptab-6-P0_S3_L003_R1_001.fastq.gz,Read1/Ptab-6-P0_S3_L004_R1_001.fastq.gz -2 Read2/Ptab-6-P0_S3_L001_R2_001.fastq.gz,Read2/Ptab-6-P0_S3_L002_R2_001.fastq.gz,Read2/Ptab-6-P0_S3_L003_R2_001.fastq.gz,Read2/Ptab-6-P0_S3_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P0_2.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-6-P5_S4_L001_R1_001.fastq.gz,Read1/Ptab-6-P5_S4_L002_R1_001.fastq.gz,Read1/Ptab-6-P5_S4_L003_R1_001.fastq.gz,Read1/Ptab-6-P5_S4_L004_R1_001.fastq.gz -2 Read2/Ptab-6-P5_S4_L001_R2_001.fastq.gz,Read2/Ptab-6-P5_S4_L002_R2_001.fastq.gz,Read2/Ptab-6-P5_S4_L003_R2_001.fastq.gz,Read2/Ptab-6-P5_S4_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P5_2.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-7-P0_S5_L001_R1_001.fastq.gz,Read1/Ptab-7-P0_S5_L002_R1_001.fastq.gz,Read1/Ptab-7-P0_S5_L003_R1_001.fastq.gz,Read1/Ptab-7-P0_S5_L004_R1_001.fastq.gz -2 Read2/Ptab-7-P0_S5_L001_R2_001.fastq.gz,Read2/Ptab-7-P0_S5_L002_R2_001.fastq.gz,Read2/Ptab-7-P0_S5_L003_R2_001.fastq.gz,Read2/Ptab-7-P0_S5_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P0_3.sam
hisat2 --max-intronlen 10000 --dta -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/HISATindex/Nb101_HISATindex/Nb101 -1 Read1/Ptab-7-P5_S6_L001_R1_001.fastq.gz,Read1/Ptab-7-P5_S6_L002_R1_001.fastq.gz,Read1/Ptab-7-P5_S6_L003_R1_001.fastq.gz,Read1/Ptab-7-P5_S6_L004_R1_001.fastq.gz -2 Read2/Ptab-7-P5_S6_L001_R2_001.fastq.gz,Read2/Ptab-7-P5_S6_L002_R2_001.fastq.gz,Read2/Ptab-7-P5_S6_L003_R2_001.fastq.gz,Read2/Ptab-7-P5_S6_L004_R2_001.fastq.gz -S /xdisk/baltrus/morgancarter/HISAT2/Ptab_nb101/Ptab_P5_3.sam
module load samtools
cd /xdisk/baltrus/morgancarter/
for i in HISAT2/Ptab_nb101/*.sam; do samtools view -bS -@ 16 $i | samtools sort -@ 16 -o $i.sorted.bam && rm $i; done
```
## StringTie v2.2.0 (HPC) - feature assignment for N. benthamiana
`StringTie_nb_B728a.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=StringTie
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH -o %j_ST_B728a_nb101.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
cd /xdisk/baltrus/morgancarter
for i in HISAT2/B728a_nb101/*.sorted.bam; do /groups/baltrus/morgan/bin/stringtie-2.2.0/stringtie $i -p 4 -G /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/Niben101.gff -o stringtie/B728a_nb101/$i.ST1.gtf -e; done
```
`StringTie_nb_Ptab.slurm`
```bash
#!/bin/bash
#SBATCH --job-name=StringTie
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH -o %j_ST_Ptab_nb101.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
cd /xdisk/baltrus/morgancarter
for i in HISAT2/Ptab_nb101/*.sorted.bam; do /groups/baltrus/morgan/bin/stringtie-2.2.0/stringtie $i -p 4 -G /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/Niben101.gff -o stringtie/Ptab_nb101/$i.ST1.gtf -e; done
```
prepDE.py3 to prepare for DESeq2 - had to move the files because they were like two more directories in then I meant.
`cd /groups/baltrus/morgan`
`python3 bin/stringtie-2.2.0/prepDE.py3 -i Ptab_B728a_RNA/plant_prepDE.txt -g Ptab_B728a_RNA/prepped/Nb101_gene_count_matrix.csv -t Ptab_B728a_RNA/prepped/Nb101_transcript_count_matrix.csv` 
```
#plant_prepDE.txt
P0_B728a_1 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P0_1.sam.sorted.bam.ST1.gtf
P0_B728a_2 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P0_2.sam.sorted.bam.ST1.gtf
P0_B728a_3 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P0_3.sam.sorted.bam.ST1.gtf
P5_B728a_1 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P5_1.sam.sorted.bam.ST1.gtf
P5_B728a_2 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P5_2.sam.sorted.bam.ST1.gtf
P5_B728a_3 /xdisk/baltrus/morgancarter/stringtie/B728a_nb101/B728a_P5_3.sam.sorted.bam.ST1.gtf
P0_Ptab_1 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P0_1.sam.sorted.bam.ST1.gtf
P0_Ptab_2 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P0_2.sam.sorted.bam.ST1.gtf
P0_Ptab_3 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P0_3.sam.sorted.bam.ST1.gtf
P5_Ptab_1 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P5_1.sam.sorted.bam.ST1.gtf
P5_Ptab_2 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P5_2.sam.sorted.bam.ST1.gtf
P5_Ptab_3 /xdisk/baltrus/morgancarter/stringtie/Ptab_nb101/Ptab_P5_3.sam.sorted.bam.ST1.gtf
```
## Bowtie v2.4.1 (HPC) - alignment for bacterial reads
Using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) that is already loaded onto the HPC. Genome information was accessed from NCBI: *P. syringae* B728a (Genbank: CP000075.1) or *P. amygdali* pv. *tabaci* (Genbank: NZ_CP042804.1 and NZ_CP042805.1) .
```
module load bowtie2
bowtie2-build Genomes/Pta_all.fa Genomes/bowtie2_index/Ptab/Ptab
bowtie2-build Genomes/Pseudomonas_syringae_B728a_112.fasta  Genomes/bowtie2_index/B728a/B728a
```
```
#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH -o %j_bowtie2_Ptab.out
#SBATCH --time=03:00:00
module load bowtie2
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_1_B5_S8_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B5_1.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_1__B0_S7_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B0_1.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_2_B0_S9_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B0_2.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_3_B0_S11_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B0_3.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_2_B5_S10_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B5_2.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/Ptab/Ptab -U /groups/baltrus/morgan/Ptab_B728a_RNA/Tabaci_reads/Ptab_3_B5_S12_R1_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/Ptab_B5_3.sam
module load samtools
for i in /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/*.sam; do samtools view -bS -@ 16 $i | samtools sort -@ 16 -o $i.sorted.bam && rm $i; done
```
```
#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH -o %j_bowtie2_B728a.out
#SBATCH --time=03:00:00
module load bowtie2
cd /groups/baltrus/morgan/Ptab_B728a_RNA/B728a_reads/
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a -1 Read1/B728a-1-B0_S1_L001_R1_001.fastq.gz,Read1/B728a-1-B0_S1_L002_R1_001.fastq.gz,Read1/B728a-1-B0_S1_L003_R1_001.fastq.gz,Read1/B728a-1-B0_S1_L004_R1_001.fastq.gz -2 Read2/B728a-1-B0_S1_L001_R2_001.fastq.gz,Read2/B728a-1-B0_S1_L002_R2_001.fastq.gz,Read2/B728a-1-B0_S1_L003_R2_001.fastq.gz,Read2/B728a-1-B0_S1_L004_R2_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B0_1.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a -1 Read1/B728a-1-B5_S4_L001_R1_001.fastq.gz,Read1/B728a-1-B5_S4_L002_R1_001.fastq.gz,Read1/B728a-1-B5_S4_L003_R1_001.fastq.gz,Read1/B728a-1-B5_S4_L004_R1_001.fastq.gz -2 Read2/B728a-1-B5_S4_L001_R2_001.fastq.gz,Read2/B728a-1-B5_S4_L002_R2_001.fastq.gz,Read2/B728a-1-B5_S4_L003_R2_001.fastq.gz,Read2/B728a-1-B5_S4_L004_R2_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B6_1.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a  -1 Read1/B728a-2-B0_S5_L001_R1_001.fastq.gz,Read1/B728a-2-B0_S5_L002_R1_001.fastq.gz,Read1/B728a-2-B0_S5_L003_R1_001.fastq.gz,Read1/B728a-2-B0_S5_L004_R1_001.fastq.gz -2 Read2/B728a-2-B0_S5_L001_R2_001.fastq.gz,Read2/B728a-2-B0_S5_L002_R2_001.fastq.gz,Read2/B728a-2-B0_S5_L003_R2_001.fastq.gz,Read2/B728a-2-B0_S5_L004_R2_001.fastq.gz  -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B0_2.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a -1 Read1/B728a-2-B5_S8_L001_R1_001.fastq.gz,Read1/B728a-2-B5_S8_L002_R1_001.fastq.gz,Read1/B728a-2-B5_S8_L003_R1_001.fastq.gz,Read1/B728a-2-B5_S8_L004_R1_001.fastq.gz -2 Read2/B728a-2-B5_S8_L001_R2_001.fastq.gz,Read2/B728a-2-B5_S8_L002_R2_001.fastq.gz,Read2/B728a-2-B5_S8_L003_R2_001.fastq.gz,Read2/B728a-2-B5_S8_L004_R2_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B5_2.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a -1 Read1/B728a-3-B0_S9_L001_R1_001.fastq.gz,Read1/B728a-3-B0_S9_L002_R1_001.fastq.gz,Read1/B728a-3-B0_S9_L003_R1_001.fastq.gz,Read1/B728a-3-B0_S9_L004_R1_001.fastq.gz -2 Read2/B728a-3-B0_S9_L001_R2_001.fastq.gz,Read2/B728a-3-B0_S9_L002_R2_001.fastq.gz,Read2/B728a-3-B0_S9_L003_R2_001.fastq.gz,Read2/B728a-3-B0_S9_L004_R2_001.fastq.gz -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B0_3.sam
bowtie2 -p 16 -x /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/bowtie2_index/B728a/B728a -1 Read1/B728a-3-B5_S12_L001_R1_001.fastq.gz,Read1/B728a-3-B5_S12_L002_R1_001.fastq.gz,Read1/B728a-3-B5_S12_L003_R1_001.fastq.gz,Read1/B728a-3-B5_S12_L004_R1_001.fastq.gz -2 Read2/B728a-3-B5_S12_L001_R2_001.fastq.gz,Read2/B728a-3-B5_S12_L002_R2_001.fastq.gz,Read2/B728a-3-B5_S12_L003_R2_001.fastq.gz,Read2/B728a-3-B5_S12_L004_R2_001.fastq.gz  -S /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/B728a_B5_3.sam
module load samtools
for i in /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/*.sam; do samtools view -bS -@ 16 $i | samtools sort -@ 16 -o $i.sorted.bam && rm $i; done
```
## featureCounts (HPC) - feature assignment for bacterial reads
FeatureCounts is part of [the subread package](http://subread.sourceforge.net/), v2.0.1 used here.
```
#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH -o %j_featureCounts_Ptab.out
#SBATCH --time=03:00:00
export PATH=$PATH:/groups/baltrus/morgan/bin/subread-2.0.1-source/bin
featureCounts \
  -T 8 \
  -t CDS \
  -g Parent \
  -M \
  - O \
  --fraction \
  -a /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/Ptab_ncbi_mod2.gff3 \
  -o /groups/baltrus/morgan/Ptab_B728a_RNA/featureCounts/featureCounts_mapped_results_Ptab.txt \
  /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/Ptab/*.sorted.bam
```
```
#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem=24gb
#SBATCH -o %j_featureCounts_Ptab.out
#SBATCH --time=03:00:00
export PATH=$PATH:/groups/baltrus/morgan/bin/subread-2.0.1-source/bin
 featureCounts \
  -p \
  -T 8 \
  -t CDS \
  -g gene_id \
  -M \
  -O \
  --fraction \
  -a /groups/baltrus/morgan/Ptab_B728a_RNA/Genomes/Pseudomonas_syringae_B728a_112.gtf \
  -o /groups/baltrus/morgan/Ptab_B728a_RNA/featureCounts/featureCounts_mapped_results_B728a.txt \
  /groups/baltrus/morgan/Ptab_B728a_RNA/bowtie2/B728a/*.sorted.bam
```
# OMA Standalone - for determining bacterial orthologs
To generate a list of genes that are orthologous between Pseudomonas syringae tabaci and Pseudomonas syringae B728a using [OMA](https://omabrowser.org/standalone/)

## Installing
Following the directions for setup from the OMA site.
```
cd /groups/baltrus/morgan/bin
wget -O oma.tgz https://omabrowser.org/standalone/OMA.2.5.0.tgz
tar xvzf oma.tgz
cd OMA.2.5.0
./install.sh /groups/baltrus/morgan/bin
cd ..
export PATH=$PATH:/groups/baltrus/morgan/bin/OMA/bin
```
## Creating Directory
```
mkdir OMA_dir
mkdir OMA_dir/DB
cp Ptab_B728a_RNA/Genomes/Psyringae_B728a_proteins.fa OMA_dir/DB/
cp Ptab_B728a_RNA/Genomes/Ptab_protein_all.fa OMA_dir/DB/
cd OMA_dir
OMA -p
OMA -c
```
The parameters file that the last command generates is how to actually change the [parameters](https://omabrowser.org/standalone/#parameters) for running OMA. You have to edit that file. I edited that and turned on `WriteOutput_PairwiseOrthologs` and added `'none'` to the `OutgroupSpecies` because that kept throwing an error.

## Running OMA
```
#!/bin/bash
#SBATCH --job-name=OMA
#SBATCH --account=baltrus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=morgancarter@email.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH -o %j_OMA.out
#SBATCH --mem=64gb
#SBATCH --time=02:00:00
cd /groups/baltrus/morgan/OMA_dir
export PATH=$PATH:/groups/baltrus/morgan/bin/OMA/bin
OMA -n 16 
```
 






> Written with [StackEdit](https://stackedit.io/).
