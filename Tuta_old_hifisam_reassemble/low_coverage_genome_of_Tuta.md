# Note for re-assemble the genome for *Tuta absoluta*
## Genome Size and Heterozygosity

- First check the assembled size for this species because it was reported as 906 Mbp  in the [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_004799115.1), and 677Mbp reported in [Tabuloc's paper](https://doi.org/10.1007/s10340-019-01116-6).

```
cat GCA_004799115.1_ASM479911v1_genomic.fna | grep -v ">" | grep [ATCGatcg] | wc -c
# 832,555,117 bps

cat GCA_004799115.1_ASM479911v1_genomic.fna | grep -v ">" | grep [ATCGatcgN] | wc -c
# 917,896,879 bps including gaps

# it is closer to the estimate on NCBI website
```

- The results of GenomeScope with different kmer size showed similar things that:
- the genome size is about 600-1000 Mbp
- the coverage is about 9X
- the heterozygosity is quite high, about 2.5% in the genome
- check the results with the links below: [K=21](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=mqt4ua4HbqnhJS5gTqjY),
[K=27](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=DVEE4sJ1cXhBHsVGw61n), and 
[K=31](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=PSyvK64NJEq2I2s1X8EF)

<br /> 

## Genome Reassembling
- A recent published paper by [Leonard et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9156671/) shows that the low read coverage (<15X) could result in high duplication rate is busco after assembling the genome with hifisam.
- This is kind of true in our case, as the coverage here is ~9X and the duplication rate is about 72%. 
- I noticed that the parameter for purging haplotig () in the script of hifisam was 2 


```
###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=hifiasm_reassemble_Tuta
#SBATCH -o hifiasm_reassemble_Tuta.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 32
#SBATCH --mem-per-cpu=10gb
#SBATCH -t 30:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara-b

module load ufrc
module load hifiasm

# to optimize the power of purging, try "-l 3"
hifiasm -o /blue/kawahara/yimingweng/Tuta_hifisam_resamble/Tuta_hifisam_resamble.asm \
        -l 3 \
        --hg-size 500m \
        --min-hist-cnt 3 \
        -t 32 \
        /blue/kawahara/shashankp/tutaabsoluta_NS2707/Tuta_genome_HiFi_assembly/m64219e_220604_102704.hifi_reads.fastq.gz

########################################################################
```

<br />


## BUSCO of the reassembled genome
- To access the quality of the reassembled genome, run busco on this genome
```
[yimingweng@login5 Tuta_reassemble]$ pwd
/blue/kawahara/yimingweng/Tuta_hifisam_reassemble/assembly1/busco_check/Tuta_reassemble

[yimingweng@login5 Tuta_reassemble]$ cat short_summary.specific.endopterygota_odb10.Tuta_reassemble.txt
# BUSCO version is: 5.3.0
# The lineage dataset is: endopterygota_odb10 (Creation date: 2020-09-10, number of genomes: 56, number of BUSCOs: 2124)
# Summarized benchmarking in BUSCO notation for file /blue/kawahara/yimingweng/Tuta_hifisam_reassemble/Tuta_hifisam_reassemble_ctg.fasta
# BUSCO was run in mode: genome
# Gene predictor used: augustus

        ***** Results: *****

        C:84.3%[S:73.3%,D:11.0%],F:1.2%,M:14.5%,n:2124
        1790    Complete BUSCOs (C)
        1557    Complete and single-copy BUSCOs (S)
        233     Complete and duplicated BUSCOs (D)
        25      Fragmented BUSCOs (F)
        309     Missing BUSCOs (M)
        2124    Total BUSCO groups searched

Dependencies and versions:
        hmmsearch: 3.1
        makeblastdb: 2.12.0+
        tblastn: 2.12.0+
        augustus: 3.4.0
        gff2gbSmallDNA.pl: None
        new_species.pl: None
        etraining: None
```

- the results is not good, but kind of acceptable, as we know that this is low (~9X) coverage genome. 
- next thing I would like to try is to rescure the data and see if I can manually make this genome better.

## Use "Purge Haplotigs pipeline" to improve the busco duplication rate
- Use [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) to clean the genome up

```
# check the dependents
module spider minimap2
# minimap2: minimap2/2.24 
module spider BEDTools
# module spider bedtools/2.30.0
module spider SAMTools
# module spider samtools/1.15
```


```
###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=Tuta_purge_haplotigs
#SBATCH -o Tuta_purge_haplotigs.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 8
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 01:00:00
#SBATCH --account=kawahara

module load minimap2
module load BEDTools
module load samtools/1.15

## step 1: map the sequence reads to the assembled genome

########################################################################
```


## Annotate the contig coverage
- So if the contigs are falsely identified from the true haplotigs, I    