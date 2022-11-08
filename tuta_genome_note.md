# Experiment note for Tuta absoluta genome assembly 
## 
The note is to record the steps we did for assembling a genome for *Tuta absoluta*. We started from a relatively small hifi read data which has only ~2.2 million reads, covering about 9-11X of the genome, assuming the genome size is around 500-600Mbp. Although we have second sequence data coming later, it was even smaller, about 0.5 million reads so it didn't help improving the assembly. Because the main issue when assembling a genome with low coverage hifi long reads is the high duplication from the duplicated haplotigs, purging the assembly without losing the read contigs is the goal. So here this note we started from purging the genome assembly.

## **10/01/2022**
I merged the two read sets, one is 9X and the other is 1-2X. I tried to get the assembly from the pooled reads. Once I got the assembly, I used haplotig purging pipeline to remove the duplicated contigs. Because multiple runs of this pipeline can potentially improve the result, I tried to run this pipeline 3 times.
1. first run 
```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_merge_purge_step1
#SBATCH -o tuta_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/hifiasm_default/Tuta_merge_default.asm.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/Tuta_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/aligned.bam  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/hifiasm_default/Tuta_merge_default.asm.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge_cutoff
#SBATCH -o tuta_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o tuta_coverage_stats.csv \
-j 80 \
-s 80
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge_haplotigs
#SBATCH -o tuta_purge_haplotigs.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/hifiasm_default/Tuta_merge_default.asm.fasta  \
-c /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/tuta_coverage_stats.csv \
-o tuta_purge_2X_50X
########################################################################
```


2. second run 

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_merge_purge_step1
#SBATCH -o tuta_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run2/tuta_purge_step2_2X_20X.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/Tuta_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run3/aligned.bam  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run2/tuta_purge_step2_2X_20X.fasta
########################################################################
```


```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge_cutoff_run2
#SBATCH -o tuta_purge_cutoff_run2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run3/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o tuta_coverage_stats_run2.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run2/tuta_purge_step2_2X_20X.fasta  \
-c /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run3/tuta_coverage_stats_run3.csv \
-o tuta_purge_step2_2X_20X
########################################################################
```
3. third run

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_merge_purge_step1
#SBATCH -o tuta_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/tuta_purge_2X_20X.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/Tuta_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run2/aligned.bam  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/tuta_purge_2X_20X.fasta
########################################################################
```


```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge_cutoff_run3
#SBATCH -o tuta_purge_cutoff_run3.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run3/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o tuta_coverage_stats_run3.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run2/tuta_purge_step2_2X_20X.fasta  \
-c /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge/purge/run3/tuta_coverage_stats_run3.csv \
-o tuta_purge_step2_2X_20X
########################################################################
```
The assembly after 3 runs of the purging pipeline still contains high duplication rate, according to the busco single copy gene evaluation.



## **10/10/2022**
Now we should give up the pooled reads, and run purge pipeline for the old single sample assembly (the one with 9X coverage).
```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge_step1
#SBATCH -o tuta_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/assembly1/Tuta_hifisam_reassemble_ctg.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/assembly1/m64219e_220604_102704.hifi_reads.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/purge/aligned.bam  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/assembly1/Tuta_hifisam_reassemble_ctg.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_original_purge_cutoff
#SBATCH -o tuta_original_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/purge/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o tuta_original_coverage_stats_run3.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/assembly1/Tuta_hifisam_reassemble_ctg.fasta  \
-c /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble/purge/tuta_original_coverage_stats_run3.csv \
-o tuta_purge_step2_2X_20X
########################################################################
```


## **10/19/2022**  
Since it's been long time without updating the progress here, I'll just quickly summarize the recent findings here.
1. The Tuta sequences from the 2 samples are neither good. The first sample has ~9-10X read coverage and the second one has only 1-2X coverage.
2. I've tried different combinations of assembling methods but non of them works perfectly. I will have to lower the standard a little bit and see what the best result I can get.
<img src="https://drive.google.com/uc?export=view&id=19HOXfstmXAP083iLs_k8vkzdXNH9PeIG">

## **10/20/2022**
**\# bwa mem**  
**\# purging**  
Because there was no perfect results from the pooled read assemblies, I would like to go back to the 9X read coverage data. At least I don't need to worry about the heterozygosity issue. Shashank used to assembled the first version of the genome with the 9X dataset and the purging setting was set to 2 (`-l 2`), and it resulted in 97% overall BUSCO but with 79% duplication rate. Then I tried rerunning hifiasm with `-l 3` to increase the power to purge the haplotigs but unfortunately I got overall BUSCO only 84%, even though the duplication rate decreased to 11%, and 3% after further purging.  
One thing I would like to try is to mapped the published short read data by [Tabuloc et al, 2019](https://link.springer.com/article/10.1007/s10340-019-01116-6) which has ~72X short read coverage, to the first version genome (the one with 97% BUSCO but high duplication rate). The reason I'd like to try this is because original step we purge the haplotigs is to map the reads we used to assemble the genome and see the coverage distribution and similarity of the contigs to determine the false splitting of haplotigs by the heterozygosity. However, since we only have 9X read depth so it is likely we don't have power to distinguish the false contigs from the alternative haplotype contigs. So let's try mapping the deep sequence reads to the genome and see if we can make the coverage cutoffs more clear. 
1. download the read files
```
[yimingweng@login2 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads

sbatch  fastq_dump.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=test_fastq_dump    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yimingweng@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=fastq_dump_test.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.9
fastq-dump --split-files --gzip  SRR8676205
########################################################################
```

2. run bwa to map the reads to the genome
```
[yimingweng@login2 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads

sbatch tuta_bwa.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=tuta_bwa
#SBATCH --output=tuta_bwa.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=180:00:00
#SBATCH -c 64

module load bwa/0.7.17
module loead samtools/1.15

bwa index /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/T_absoluta_hifiasm_06_14_2022.asm.bp.p_ctg.fa

bwa mem -t 64 /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/T_absoluta_hifiasm_06_14_2022.asm.bp.p_ctg.fa /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/SRR8676205_2.fastq.gz /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/SRR8676205_3.fastq.gz > tuta_short_read_aln.sam

samtools view -S -b tuta_short_read_aln.sam > tuta_short_read_aln.bam
samtools sort -m 1G -o tuta_short_read_aln.bam -T tmp.ali
########################################################################
```

## **10/31/2022** 
**\# purging** 

Once we have the mapped bam file, use it to run the haplotig purging again
```
[yimingweng@login5 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads

sbatch tuta_purge2_step2.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_purge2_step2
#SBATCH -o tuta_purge2_step2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_short_read_aln_sorted.bam \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/tuta_purge.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=tuta_original_purge_cutoff
#SBATCH -o tuta_original_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_short_read_aln_sorted.bam.gencov  \
-l 30  \
-m 85  \
-h 190  \
-o tuta_run2_stats.csv \
-j 95 \
-s 95

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/tuta_purge.fasta  \
-c /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_run2_stats.csv \
-o tuta_purge_run2
########################################################################
Result: C:96.2%[S:82.5%,D:13.7%],F:0.5%,M:3.3%,n:5286
```


## **11/02/2022** 
**\# best assembly model**  
**\# blobplot**  

By considering the BUSCO score, the best assembly I can get is to use the 9X genomic data, and run hifiasm with purging aggressiveness to be `-l 2`, and let the haplotig-purging pipeline to remove the duplications. So here is the final assembly for annotation:`/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta`

1. Check the assembly statistics:
```
[yimingweng@login5 purge]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge

module load python3

python /blue/kawahara/yimingweng/universal_scripts/assemblystats.py /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta
{
  "Contig Stats": {
    "L10": 12,
    "L20": 29,
    "L30": 51,
    "L40": 79,
    "L50": 115,
    "N10": 4435971,
    "N20": 3284839,
    "N30": 2611967,
    "N40": 2067342,
    "N50": 1614219,
    "gc_content": 38.454626166301814,
    "longest": 7993428,
    "mean": 948696.4491279069,
    "median": 630478.0,
    "sequence_count": 688,
    "shortest": 3255,
    "total_bps": 652703157
```
-> L50: **1.61M**  
-> Assembled size" **65 Mbps**

2. run blobplot to identify the potential foreign sequences
    - run megablast to assign the best hit gene to the contig
```
[yimingweng@login5 blobplot]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/blobplot

sbatch -J tuta_genome /blue/kawahara/yimingweng/universal_scripts/megablast_nt.slurm /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta tuta_genome
```

3. Store this "final version" of the assembly in `/blue/kawahara/yimingweng/Tuta_genome_project/assemblies` and name it as "**tuta_final_assembly.fasta**".
```
yimingweng@login2 assemblies]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/assemblies

cp /blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_sho
rt_reads/tuta_purge_run2.fasta  ./

mv tuta_purge_run2.fasta tuta_final_assembly.fasta
```




```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
. blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp


sbatch -J tuta /blue/kawahara/yimingweng/universal_scripts/minimap.slurm \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta \
/blue/kawahara/shashankp/tutaabsoluta_NS2707/Tuta_fastqc/m64219e_220604_102704.hifi_reads.fasta.gz \
tuta_blobplot


sbatch -J tuta /blue/kawahara/yimingweng/universal_scripts/blobplot.slurm \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/blobplot/tuta_blobplot.bam \
/blue/kawahara/yimingweng/Tuta_genome_project/blobplot/tuta_genome.nt.mts1.hsp1.1e25.megablast.out \
tuta_blobplot

sbatch -J tuta /blue/kawahara/yimingweng/universal_scripts/repeatmodeler2.slurm \
/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta_purge_run2.fasta \
tuta_repeat
```
<img src="https://github.com/yimingweng/Tuta_genome_project/blob/main/blobplot_results/tuta_blobplot.png">


```
[yimingweng@login2 repeatmasker]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/annotation/repeatmasker

sbatch -J tuta /blue/kawahara/yimingweng/universal_scripts/repeatmakser.slurm /blue/kawahara/yimingweng/Tuta_genome_project/assemblies/tuta_final_assembly.fasta \
/blue/kawahara/yimingweng/Tuta_genome_project/annotation/repeatmodeler/tuta_repeat-families.fa \
tuta_repeatmasker

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_repeatmask_%j
#SBATCH -o %x_repeatmask_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 32

genome=${1}
modeler_out={2}
prefix=${3}

mkdir ${prefix}

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1


name=$(echo ${genome} | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
path=$(pwd)

# step 1: mask the simple repeat
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-noint \
-no_is \
-dir ${prefix} \
${genome} &> ./${prefix}/${prefix}_step1.out


# step 2: mask repeats based on existing databases
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-species Lepidoptera \
-dir ${prefix} \
${path}/${prefix}/${name}.masked &> ./${prefix}/${prefix}_step2.out


# step 3: mask genome based on the output of repeatmodeler
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-lib ${modeler_out} \
-dir ${prefix} \
${path}/${prefix}/${name}.masked.masked &> ./${prefix}/${prefix}_step3.out
```

```
[yimingweng@login6 tuta_repeatmasker]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/annotation/repeatmasker/tuta_repeatmasker

bash /blue/kawahara/yimingweng/universal_scripts/maskrate.sh tuta_final_assembly.fasta.masked.masked.masked
softmasking rate is 54.40%
hardmasking rate is 0%
```


braker2
```
[yimingweng@login6 prothint]$ pwd
/blue/kawahara/yimingweng/Tuta_genome_project/annotation/braker2/prothint

sbatch -J tuta /blue/kawahara/yimingweng/universal_scripts/braker_prothint.slurm \
/blue/kawahara/yimingweng/Tuta_genome_project/annotation/repeatmasker/tuta_repeatmasker/tuta_final_assembly.fasta.masked.masked.masked \
tuta_prothint
```
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  


### **Path change/File moving notes**
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_merge`
- 11/02/2022: rename directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_new_seq/` to be `~/Tuta_second_sequence/`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_second_sequence/genome_size`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_old_hifisam_reassemble`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/*fastq.gz`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Tuta_genome_project/Tuta_first_assembly/purge/run2_with_short_reads/tuta*.bam`