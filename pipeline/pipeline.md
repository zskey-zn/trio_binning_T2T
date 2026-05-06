# 00.trio_binning_assembly
Use hifiasm/verkko2 generate Trio binning contigs.
## hifiasm

```bash
HiFi=all_hifi.bam.fasta
ONT=all.sup.pass.fq.gz
pat_1=paternal_1.clean.fq.gz
pat_2=paternal_2.clean.fq.gz
mat_1=maternal_1.clean.fq.gz
mat_2=maternal_2.clean.fq.gz
prefix=Sus_scrofa

yak count -b37 -t16 -o pat.yak <(cat ${pat_1} ${pat_2}) <(cat ${pat_1} ${pat_2})
yak count -b37 -t16 -o mat.yak <(cat ${mat_1} ${mat_2}) <(cat ${mat_1} ${mat_2})

hifiasm -o ${prefix} -t 50 -1 ./pat.yak -2 ./mat.yak --ul ${ONT} ${HiFi} 2> out.log
awk '/^S/{print ">"$2;print $3}' ${prefix}.dip.hap1.p_ctg.gfa > ${prefix}.dip.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${prefix}.dip.hap2.p_ctg.gfa > ${prefix}.dip.hap2.p_ctg.fa
```

## verkko

```bash
# generate compress kmer database
meryl count k=30 threads=24 compress output Mother.meryl maternal_*.clean.fq.gz
meryl count k=30 threads=24 compress output father.meryl paternal_*.clean.fq.gz
HiFi=all_hifi.bam.fasta
ONT=all.sup.pass.fq.gz
# run verkko2 
verkko -d asm \
  --hifi $HiFi \
  --nano $ONT \
  --threads 30 --slurm \
  --hap-kmers father.hapmer.meryl \
              Mother.hapmer.meryl \
              trio
```

## 01.reads_split
Use `yak triobin` to genotype HiFi/ONT reads.
<img width="563" height="191" alt="image" src="https://github.com/user-attachments/assets/f3a417f0-2e68-4f2f-9165-75c16bacc249" />
Parameter reference：https://github.com/lh3/yak/issues/1
```bash

```
