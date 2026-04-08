cat ${sample}_partial_insertion_200kb.txt |sed 's/\[//g' |sed 's/\]//g' |sed "s/'//g" |sed 's/, /\t/g' > temp && mv temp  ${sample}_partial_insertion_200kb.txt
cat ${sample}_spanning_insertion.txt |sed 's/\[//g' |sed 's/\]//g' |sed "s/'//g" |sed 's/, /\t/g' > temp && mv temp ${sample}_spanning_insertion.txt

awk '{OFS="\t"}{print $1,$2,$3,$6"."$4"."$5"."$7,$8,"p"}' ${sample}_partial_insertion_200kb.txt > ${sample}_partial_insertion_200kb.tmp
awk '{OFS="\t"}{print $2,int(int(($3+$4)/2)/200)*200,int(int(($3+$4)/2)/200)*200+200,int(($3+$4)/2)"."$5"."$6"."$7"_"$8,$1,"s"}' ${sample}_spanning_insertion.txt > ${sample}_spanning_insertion_200kb.tmp
source /sibcb/program/install/python-3.9/profile
python ~/Data/WangLu/20240411_E13.5-Cramp1-KO-DNA-nanopore/20240512_NewSum.py ${sample}_partial_insertion_200kb.tmp ${sample}_spanning_insertion_200kb.tmp ${sample}_Summary.bed
rm ${sample}_partial_insertion_200kb.tmp ${sample}_spanning_insertion_200kb.tmp
awk '{if ($5+$6 <=2) print $0}' ${sample}_Summary.bed > ${sample}_SupportReads2.bed
