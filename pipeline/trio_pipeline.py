# -*- coding: utf-8 -*-
"""
Created on Fri May  8 11:04:47 CST 2026

@author: zhengshang

"""
import re
import sys
import os
data_from=sys.argv[1].strip()
work_dir=os.getcwd()
##conf reading
f=open(data_from,'r')
data=f.readlines()
f.close()
data=[i[:-1] for i in data if i[0]!="#" and i!="\n"]
s_a=[i.split('=')[0].strip() for i in data]
s_b=[i.split('=')[1].strip() for i in data]
data_dict=dict(zip(s_a,s_b))

os.chdir(work_dir)
os.system('mkdir 00.data 01.assembly_all 02.assembly_separately')
os.chdir('./00.data')
os.system('mkdir yak meryl')
os.chdir('./yak')
f=open('step1.generate_parental_yak.sh','w')
f.write("""#!/bin/bash
mat_1={}
mat_2={}
pat_1={}
pat_2={}""".format(data_dict['mat_1'],data_dict['mat_2'],data_dict['pat_1'],data_dict['pat_2'])+"""
/public/frasergen/PUB/software/yak/yak count -b37 -t16 -o pat.yak <(cat ${pat_1} ${pat_2}) <(cat ${pat_1} ${pat_2})
/public/frasergen/PUB/software/yak/yak count -b37 -t16 -o mat.yak <(cat ${mat_1} ${mat_2}) <(cat ${mat_1} ${mat_2})
touch step1.generate_parental_yak.sh.finish """)
f.close()
f=open('step2.hifi_genotype.sh','w')
f.write("""#!/bin/bash
module load seqkit
hifireads={}""".format(data_dict['hifireads'])+"""
/public/frasergen/PUB/software/yak/yak triobin -t 16  pat.yak mat.yak $hifireads > hifi_triobin.txt
awk '$3>=21&&$4<=2&&$2=="p"{print $1}' hifi_triobin.txt > hifi_paternal.txt
awk '$4>=21&&$3<=2&&$2=="m"{print $1}' hifi_triobin.txt > hifi_maternal.txt
seqkit grep -f hifi_paternal.txt $hifireads >hifi_pat.fa
seqkit grep -f hifi_maternal.txt $hifireads >hifi_mat.fa 
touch step2.hifi_genotype.sh.finish""")
f.close()
f=open('step3.ont_genotype.sh','w')
f.write("""#!/bin/bash
module load seqkit
ontreads={}""".format(data_dict['ontreads'])+"""
/public/frasergen/PUB/software/yak/yak triobin -t 16  pat.yak mat.yak $ontreads > ont_triobin.txt
awk '$3>=21&&$4<=2&&$2=="p"{print $1}' ont_triobin.txt > ont_paternal.txt
awk '$4>=21&&$3<=2&&$2=="m"{print $1}' ont_triobin.txt > ont_maternal.txt
seqkit grep -f ont_paternal.txt $ontreads >ont_pat.fq
seqkit grep -f ont_maternal.txt $ontreads >ont_mat.fq
touch step3.ont_genotype.sh.finish """)
f.close()

os.chdir('../meryl')
f=open("step1_meryl_father.sh","w")
f.write("""#!/bin/bash
export PATH=/public/frasergen/PUB/software/meryl-1.4/meryl-1.4.1/bin:$PATH
meryl count k=30 threads=24 compress output father.meryl {} {}
touch step1_meryl_father.sh.finish""".format(data_dict['pat_1'],data_dict['pat_2']))
f.close()
f=open("step2_meryl_mother.sh","w")
f.write("""#!/bin/bash
export PATH=/public/frasergen/PUB/software/meryl-1.4/meryl-1.4.1/bin:$PATH
meryl count k=30 threads=24 compress output mother.meryl {} {}
touch step2_meryl_mother.sh.finish""".format(data_dict['mat_1'],data_dict['mat_2']))
f.close()
f=open("step3_meryl_child.sh","w")
f.write("""#!/bin/bash
export PATH=/public/frasergen/PUB/software/meryl-1.4/meryl-1.4.1/bin:$PATH
meryl count k=30 threads=24 compress output child.meryl {} {}
touch step3_meryl_child.sh.finish""".format(data_dict['R1'],data_dict['R2']))
f.close()
f=open("step4.hapmers.sh","w")
f.write("""#!/bin/bash
export PATH=/public/frasergen/PUB/software/meryl-1.4/meryl-1.4.1/bin:$PATH
export PATH=/public/frasergen/PUB/software/merqury-1.4:$PATH
/public/frasergen/PUB/software/merqury-1.4/trio/hapmers.sh  Mother.meryl father.meryl child.meryl
touch step4.hapmers.sh.finish""")
f.close()

os.chdir('../../01.assembly_all')
os.system('mkdir 00.hifiasm_hic 01.hifiasm_trio 02.verkko2_trio')
os.chdir('./00.hifiasm_hic')
f=open('step00.hifiasm_hic.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}
HiFi={}
ONT={}
hic1={}
hic2={}
telomere_motif={} """.format(data_dict['sp'],data_dict['hifireads'],data_dict['ontreads'],data_dict['hic1'],data_dict['hic2'],data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} --telo-m ${telomere_motif} -t 50 --ul ${ONT} --h1 ${hic1} --h2 ${hic2} ${HiFi} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.hic.hap1.p_ctg.gfa > ${prefix}.hic.hap1.p_ctg.fa && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.hic.hap2.p_ctg.gfa > ${prefix}.hic.hap2.p_ctg.fa && \\
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.hic.hap*.p_ctg.fa >> N50 && \\

echo end at `date`
touch step00.hifiasm_hic.sh.finish """)
f.close()
os.chdir('../01.hifiasm_trio')
f=open('step01.hifiasm_trio.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}
HiFi={}
ONT={}
pat={}/00.data/yak/pat.yak
mat={}/00.data/yak/mat.yak
telomere_motif={} """.format(data_dict['sp'],data_dict['hifireads'],data_dict['ontreads'],work_dir,work_dir,data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} -1 ${pat} -2 ${mat} --telo-m ${telomere_motif} -t 50 --ul ${ONT} ${HiFi} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.dip.hap1.p_ctg.gfa > ${prefix}.dip.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${prefix}.dip.hap2.p_ctg.gfa > ${prefix}.dip.hap2.p_ctg.fa
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.dip.hap*.p_ctg.fa > N50

echo end at `date`
touch step01.hifiasm_trio.sh.finish
""")
f.close()
os.chdir('../02.verkko2_trio')
f=open('step02.verkko2_trio.sh','w')
f.write("""source  /public/frasergen/PUB/software/conda/miniconda3/bin/activate
conda activate /public/frasergen/PUB/software/conda/miniconda3/envs/verkko_20241204
export PATH=/public/frasergen/PUB/software/conda/miniconda3/envs/verkko_20241204/lib/verkko/bin:$PATH
export LD_LIBRARY_PATH=/public/frasergen/PUB/software/conda/miniconda3/envs/verkko_20241204/lib:$LD_LIBRARY_PATH

HiFi={}
ONT={}
pat={}/00.data/meryl/father.hapmer.meryl
mat={}/00.data/meryl/Mother.hapmer.meryl
telomere_motif={} """.format(data_dict['hifireads'],data_dict['ontreads'],work_dir,work_dir,data_dict['telomere_motif'])+"""

verkko -d asm \\
  --hifi $HiFi \\
  --nano $ONT \\
  --threads 30 --telomere-motif ${telomere_motif} --slurm \\
  --hap-kmers ${pat} \\
              ${mat} \\
              trio

touch step02.verkko2_trio.sh.finish
""")
f.close()

os.chdir('../../02.assembly_separately')
os.system('mkdir 01.ul 02.ont')
os.chdir('./01.ul')
os.system('mkdir 01.pat 02.mat')
os.chdir('./01.pat')
f=open('hifiasm_ul_pat.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}_ul_pat
HiFi={}/00.data/yak/hifi_pat.fa
ONT={}/00.data/yak/ont_pat.fq
telomere_motif={} """.format(data_dict['sp'],work_dir,work_dir,data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} --telo-m ${telomere_motif} -t 26 --ul ${ONT} ${HiFi} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.bp.p_ctg.gfa > ${prefix}.bp.p_ctg.fa
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.bp.p_ctg.fa >> N50 && \\

echo end at `date`
touch hifiasm_ul_pat.sh.finish
""")
f.close()
os.chdir('../02.mat')
f=open('hifiasm_ul_mat.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}_ul_mat
HiFi={}/00.data/yak/hifi_mat.fa
ONT={}/00.data/yak/ont_mat.fq
telomere_motif={} """.format(data_dict['sp'],work_dir,work_dir,data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} --telo-m ${telomere_motif} -t 26 --ul ${ONT} ${HiFi} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.bp.p_ctg.gfa > ${prefix}.bp.p_ctg.fa
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.bp.p_ctg.fa >> N50 && \\

echo end at `date`
touch hifiasm_ul_mat.sh.finish
""")
f.close()
os.chdir('../../02.ont')
os.system('mkdir 01.pat 02.mat')
os.chdir('./01.pat')
f=open('hifiasm_ont_pat.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}_ont_pat
ONT={}/00.data/yak/ont_pat.fq
telomere_motif={} """.format(data_dict['sp'],work_dir,data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} --telo-m ${telomere_motif} -t 26 --ul ${ONT} --ont ${ONT} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.bp.p_ctg.gfa > ${prefix}.bp.p_ctg.fa
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.bp.p_ctg.fa >> N50 && \\

echo end at `date`
touch hifiasm_ont_pat.sh.finish
""")
f.close()
os.chdir('../02.mat')
f=open('hifiasm_ont_mat.sh','w')
f.write("""#!/bin/bash
echo begin at `date`
prefix={}_ont_mat
ONT={}/00.data/yak/ont_mat.fq
telomere_motif={} """.format(data_dict['sp'],work_dir,data_dict['telomere_motif'])+"""

/public/frasergen/PUB/software/hifiasm/hifiasm-0.25.0/hifiasm -o ${prefix} --telo-m ${telomere_motif} -t 26 --ul ${ONT} --ont ${ONT} 2> out.log && \\
awk '/^S/{print ">"$2;print $3}' ${prefix}.bp.p_ctg.gfa > ${prefix}.bp.p_ctg.fa
/usr/bin/java -jar /public/frasergen/PUB/software/gnx-tools/gnx-tools-master/gnx.jar ${prefix}.bp.p_ctg.fa >> N50 && \\

echo end at `date`
touch hifiasm_ont_mat.sh.finish
""")
f.close()

os.chdir(work_dir)
f=open('shell.list','w')
f.write("""{0}/00.data/meryl/step1_meryl_father.sh:100g:24
{0}/00.data/meryl/step2_meryl_mother.sh:100g:24
{0}/00.data/meryl/step3_meryl_child.sh:100g:24
{0}/00.data/yak/step1.generate_parental_yak.sh:64g:16
{0}/01.assembly_all/00.hifiasm_hic/step00.hifiasm_hic.sh:300g:50
{0}/00.data/meryl/step1_meryl_father.sh:100g:24\t{0}/00.data/meryl/step4.hapmers.sh:12g:3
{0}/00.data/meryl/step2_meryl_mother.sh:100g:24\t{0}/00.data/meryl/step4.hapmers.sh:12g:3
{0}/00.data/meryl/step3_meryl_child.sh:100g:24\t{0}/00.data/meryl/step4.hapmers.sh:12g:3
{0}/00.data/yak/step1.generate_parental_yak.sh:64g:16\t{0}/00.data/yak/step2.hifi_genotype.sh:64g:16
{0}/00.data/yak/step1.generate_parental_yak.sh:64g:16\t{0}/00.data/yak/step3.ont_genotype.sh:64g:16
{0}/00.data/yak/step1.generate_parental_yak.sh:64g:16\t{0}/01.assembly_all/01.hifiasm_trio/step01.hifiasm_trio.sh:300g:50
{0}/00.data/meryl/step4.hapmers.sh:12g:3\t{0}/01.assembly_all/02.verkko2_trio/step02.verkko2_trio.sh:120g:30
{0}/00.data/yak/step3.ont_genotype.sh:64g:16\t{0}/02.assembly_separately/02.ont/01.pat/hifiasm_ont_pat.sh:120g:26
{0}/00.data/yak/step3.ont_genotype.sh:64g:16\t{0}/02.assembly_separately/02.ont/02.mat/hifiasm_ont_mat.sh:120g:26
{0}/00.data/yak/step3.ont_genotype.sh:64g:16\t{0}/02.assembly_separately/01.ul/01.pat/hifiasm_ul_pat.sh:120g:26
{0}/00.data/yak/step3.ont_genotype.sh:64g:16\t{0}/02.assembly_separately/01.ul/02.mat/hifiasm_ul_mat.sh:120g:26
{0}/00.data/yak/step2.hifi_genotype.sh:64g:16\t{0}/02.assembly_separately/01.ul/01.pat/hifiasm_ul_pat.sh:120g:26
{0}/00.data/yak/step2.hifi_genotype.sh:64g:16\t{0}/02.assembly_separately/01.ul/02.mat/hifiasm_ul_mat.sh:120g:26
""".format(work_dir))
f.close()
f=open('sbatch.sh','w')
f.write("""#!/bin/bash
perl /public/frasergen/DNA/pipeline/02.assembly/05.report/script/Monitor_SGE_v0.4.pl -list {0}/shell.list -q xhacexclu11,xhacexclu16,xhacexclu03 -o {0}/job.log
""".format(work_dir))
f.close()



