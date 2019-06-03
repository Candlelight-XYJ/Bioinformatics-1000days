## 学习材料
+ https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/
## 1. VCF（Variant Calling Format）
VCF文件是常见的文本格式，它包含三个主要内容：
+ meta-information lines（元信息）
+ a head line（列名信息）
+ data lines（包含基因组中某个位点的信息，也包括基因型的信息）

```r
##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3

```

这个##符号后解释的信息比较多，可以自己查看手册，这里重点讲解每一列的含义
每份VCF文件有8列是确定的，每一列之间都是以tab分隔符隔开，且缺失数据的位置一般使用**`.`** 符号来替代数据

+ **CHROM-chromosome**：染色体
+ **POS-position**：参考位点的位置，基因组上起始的第一个碱基的位置为1。在每个参考序列的染色体中，位点位置按递增顺序进行数值排序。允许有多个POS相同的记录，端粒用位置0或N+1表示，其中N为对应染色体或重叠群的长度。(要求N为整数)
+ **ID-identifier**：记录的ID号，可以是SNP的rsid，如果没有id号，每个id号只在一条记录中显示。
+ **REF-reference base(s)**：参考碱基，碱基为A,C,T,G,N中的任一个（不区分大小写），一条记录中的ref列可以有多个碱基，**`此时POS指的就是碱基字符串中第一个碱基的位置`**
+ **ALT-alternate bases(s)** ：由逗号分隔的,相对于参考基因组的(变异)碱基
+ **QUAL - quality**： Phred格式(Phred_scaled)的质量值，可以理解为所call出来的变异位点的质量值。表示在该位点存在variant的可能性；`该值越高，则variant的可能性越大； `
计算方法：① Q=-10*lgP，Q表示质量值；P表示这个位点发生错误的概率。 
②Phred值Q = -10 * lg (1-p) ，p为variant存在的概率; 
通过计算公式可以看出值为10的表示错误概率为0.1，该位点为variant的概率为90%。 
同理，当Q=20时，错误率就控制在了0.01。
+ **FILTER - filter status**：变异位点的过滤记录。FILTER的PASS代表变异位点通过了过滤，是比较好的标准变异。如果没有通过过滤，就会在FILTER这一栏提示除了PASS的其他信息。如果这一栏是一个“.”的话，就说明没有进行过任何过滤。
+ **INFO**：其它信息-较多