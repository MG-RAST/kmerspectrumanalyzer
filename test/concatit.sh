cd total

ln -s ../orig/SRR036083.fastq SRR036083.fastq
ln -s ../orig/SRR035445.fastq SRR035445.csfastq
ln -s ../orig/SRR387449.fastq SRR387449.fastq
ln -s ../orig/SRR341577.fastq SRR341577.fastq
ln -s ../orig/SRR006331.fastq SRR006331.fastq
ln -s ../orig/SRR000333.fastq SRR000333.fastq 
ln -s ../orig/SRR016816.fastq SRR016816.fastq   # late
ln -s ../orig/SRR402885.fastq SRR402885.fastq   # late

cd ..

cat orig/SRR036211_?.fastq > total/SRR036211B.fastq
cat orig/SRR071425_?.fastq > total/SRR071425B.fastq
cat orig/SRR089543_?.fastq > total/SRR089543B.fastq
cat orig/SRR089544_?.fastq > total/SRR089544B.fastq
cat orig/SRR090599_?.fastq > total/SRR090599B.fastq
cat orig/SRR190843_?.fastq > total/SRR190843B.fastq
cat orig/SRR387449_?.fastq > total/SRR387449B.fastq
cat orig/SRR610299_?.fastq > total/SRR610299B.fastq
cat orig/SRR610309_?.fastq > total/SRR610309B.fastq
cat orig/SRR769599_?.fastq > total/SRR769599B.fastq
cat orig/SRR769600_?.fastq > total/SRR769600B.fastq
cat orig/SRR769601_?.fastq > total/SRR769601B.fastq
cat orig/SRR769602_?.fastq > total/SRR769602B.fastq
cat orig/SRR769603_?.fastq > total/SRR769603B.fastq
cat orig/SRR488146_?.fastq > total/SRR488146B.fastq  # late

cat orig/SRR03996[3456].fastq > total/SRR039966A.fastq 
cat orig/SRR006330.fastq orig/SRR006332.fastq > total/SRR006332A.fastq
cat orig/SRR072776.fastq orig/SRR090721.fastq > total/SRR090721A.fastq 
cat orig/SRR029227.fastq orig/SRR000333.fastq > total/SRR029227A.fastq 
cat orig/SRR059788.fastq orig/SRR059789.fastq > total/SRR059789A.fastq 
cat orig/SRR004368.fastq orig/SRR017401.fastq > total/SRR017401A.fastq
cat orig/SRR060959.fastq orig/SRR060960.fastq > total/SRR060960A.fastq
cat orig/SRR072318.fastq orig/SRR090709.fastq > total/SRR090709A.fastq
cat orig/SRR089543.fastq orig/SRR089544.fastq > total/SRR089544A.fastq
cat orig/SRR059232.fastq orig/SRR059233.fastq orig/SRR059234.fastq orig/SRR059235.fastq orig/SRR059236.fastq > total/SRR059236A.fastq
cat orig/SRR396647.fastq orig/SRR396649.fastq orig/SRR396650.fastq orig/SRR396651.fastq orig/SRR396653.fastq > total/SRR396651A.fastq 
cat orig/SRR031261.fastq orig/SRR031262.fastq orig/SRR031263.fastq orig/SRR031264.fastq orig/SRR031265.fastq orig/SRR031266.fastq > total/SRR031266A.fastq # mistake cat first time



