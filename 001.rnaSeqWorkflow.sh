#!/bin/bash
#
#本脚本使用前提，PE数据，格式如下：SRR20736630_1.fastq.gz SRR20736630_2.fastq.gz ......等
#fastp shell script
if [ ! -d cleandata  ]; then
	mkdir cleandata
elif [ ! -d html ]; then
	mkdir html
elif [ ! -d json ]; then
	mkdir json
fi

ls *.gz > list 
for i in `cat list`;do
	sampleid=${i:0:11} 			#实操过程中这里可能要改动一下
	echo ${sampleid} >> list1
	echo ${sampleid}"_1.fastq.gz" >> list2
	echo ${sampleid}"_2.fastq.gz" >> list3
done
paste list1 list2 list3 > list
sort list | uniq > list1 ; mv list1 list

awk '{print "fastp -w 30 -q 30 -i "$2" -o ./cleandata/"$2" -I "$3" -O ./cleandata/"$3" -h ./html/"$1".html -j ./json/"$1".json"}' list > fastp.sh
#bash fastp.sh
awk '{print "mkdir "$1}' list > mkdirforjson.sh
#bash mkdirforjson.sh
awk '{print "cp ./json/"$1".json ./"$1"/fastp.json"}' list > cpjson.sh #驱动multiqc的先决条件
#bash cpjson.sh

#multiqc .

awk '{print "rm -r "$1}' list > rmdirforjson.sh
#bash rmdirforjson.sh

rm list2 list3 cpjson.sh mkdirforjson.sh fastp.sh rmdirforjson.sh

#================================================================================================

mkdir raw_data
# mv *.gz raw_data
#cp ./cleandata/*.gz .
nowpath=`pwd`
refpath="/share/liuqi/songrunxian/0.liuqi/experiment_010_medicinekeyan/ref/HADb"
fasta="hg38.fa"
gtf="human_HADb.gtf"
cd ${refpath}
echo "STAR --runMode genomeGenerate --runThreadN 20 --genomeDir ./index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf}" >> starindex.sh
#bash starindex.sh
cd ${nowpath}
echo "for i in \`cut -f 1 list\`;do
	STAR --genomeDir ${refpath}/index --readFilesCommand zcat --runThreadN 15 --outSAMtype BAM SortedByCoordinate --outSAMattributes MD NH --readFilesIn \${i}_1.fastq.gz \${i}_2.fastq.gz --outFileNamePrefix \${i}
done" > star.sh
#bash star.sh
rm star.sh
rename Aligned.sortedByCoord.out.bam .bam *
ls *.bam > bamlist
#featureCounts -T 10 -a ${refpath}/${gtf} -o allsample_read.count -p -B -C -f -t exon -g gene_id `cut -f 1 bamlist`
rm bamlist

awk '{if(($7!=0)||($8!=0)||($9!=0)||($10!=0)||($11!=0)||($12!=0)) print $0}' allsample_read.count > simple_allsample_read.count

cut -f 1,7-100 simple_allsample_read.count > sim_allsample_read.count
cut -f 1,6 simple_allsample_read.count > gene.lengths.txt

echo "library(tidyverse)
library(readr)
file <- read_delim(\"sim_allsample_read.count\", delim=\"\\t\", escape_double=FALSE, comment=\"#\", trim_ws=TRUE)
sum_file <- file %>%
  group_by(Geneid) %>%
  summarise(SRR20736630=sum(SRR20736630.bam),
            SRR20736631=sum(SRR20736631.bam),
            SRR20736632=sum(SRR20736632.bam),
            SRR20736633=sum(SRR20736633.bam),
            SRR20736634=sum(SRR20736634.bam),
            SRR20736635=sum(SRR20736635.bam)
            )
counts <- sum_file
write.table(counts, file = \"count.txt\", sep = \"\\t\", row.names = FALSE)
counts <- read.table(\"count.txt\", header=TRUE, row.names=1, stringsAsFactors=FALSE)

gene_lengths <- read_delim(\"gene.lengths.txt\", delim = \"\\t\", escape_double = FALSE,comment = \"#\", trim_ws = TRUE)
genelen <- gene_lengths%>%
  group_by(Geneid) %>%
  summarise(length=sum(Length))
write.table(genelen, file = \"gl.txt\", sep = \"\\t\", row.names = FALSE)
gene_lengths <- read.table(\"gl.txt\", header=TRUE, row.names=1, stringsAsFactors=FALSE)

gene_lengths <- gene_lengths[rownames(counts),] #按照counts中行名，排列成向量列表
total_reads <- colSums(counts)
fpkm <- t(t(counts) / gene_lengths) / (total_reads / 1e6)
#write.table()函数的quote参数用于指定是否在输出文件中使用引号。如果将其设置为TRUE，则所有数据都将用引号括起来。如果将其设置为FALSE，则不会使用引号
write.table(fpkm, file=\"fpkm.txt\", sep=\"\\t\", quote=FALSE)" >> fpkm.r
#Rscript fpkm.r
mv fpkm.txt allsample_fpkm.txt
rm fpkm.r gl.txt count.txt gene.lengths.txt simple_allsample_read.count sim_allsample_read.count

echo "library(GEOquery)
gset <- getGEO(\"GSE210249\", getGPL=FALSE)
gset <- gset[[1]]
sample_info <- pData(gset)
write.table(sample_info, file=\"sample_info.txt\", sep=\"\\t\", quote=FALSE)" > mksample_info.r
#Rscript mksample_info.r
rm mksample_info.r
#此时需要修改一下样本表
echo "library(readr)
sample_info <- read_delim(\"sample_info.txt\", 
    delim = \"\t\", escape_double = FALSE, 
    trim_ws = TRUE)
library(tidyverse)
spin <- sample_info
s <- spin %>% 
  select(geo_accession, title) %>%
  mutate(group = str_split(title, '-', simplify = T)[,2]) %>%
  mutate(group_num = case_when(
    group == 'SP' ~ 1,
    group == 'AD' ~ 2 )) %>%
  mutate(test1 = round(runif(6,1,100), digits = 2),
         test2 = round(runif(6,1,100), digits = 2),
         test3 = round(runif(6,1,100), digits = 2),
         test4 = round(runif(6,1,100), digits = 2 )) %>% 
  column_to_rownames(var = 'geo_accession') 

group_by(s, group) %>%
  summarise(count = n())

write.table(s, file=\"sample_information.txt\", sep=\"\\t\", quote=FALSE)
    " > mksaple_info.r
#Rscript mksaple_info.r
rm mksaple_info.r sample_info.txt
#手动下载基因注释表，在GEO浏览器中
mkdir mindata
cp allsample_fpkm.txt mindata
cp sample_information.txt mindata
cp Human.GRCh38.p13.annot.tsv mindata
tar -zvcf mindata.tar.gz mindata
#================================================================================================



