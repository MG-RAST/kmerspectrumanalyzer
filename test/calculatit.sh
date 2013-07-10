

cd total
for i in *.fastq; do echo $i; countkmer21.sh $i; done 

mkdir fits
cd fits
ln ../total/*.21 .


for i in *.21 
do 
if [[ ! -e $i.fit.png ]] 
then
kmerspectrumanalyzer.py -q -n 30 $i 
fi 
done 
