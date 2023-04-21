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
bash mkdirforjson.sh
awk '{print "cp ./json/"$1".json ./"$1"/fastp.json"}' list > cpjson.sh
bash cpjson.sh

multiqc .

awk '{print "rm -r "$1}' list > rmdirforjson.sh
bash rmdirforjson.sh

rm list2 list3 cpjson.sh mkdirforjson.sh fastp.sh rmdirforjson.sh
