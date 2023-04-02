#This isi to check if any mmpbsa runs have missed or need to be repeated

cd /home/chemical/phd/chz198152/scratch/covid/kb1-chC1/mmpbsa/
<< 'END'
sort -u mmpbsa.txt > uniq_mmpbsa.txt
sed -i '/^$/d' uniq_mmpbsa.txt

#1. Find completely  missed runs
e=`wc -l names.txt | awk '{print $1}'`
for ((j=1; j<=$((e));  j++))
do

l=`sed "${j}q;d" names.txt`

if grep -Fq "$l" uniq_mmpbsa.txt
then
	echo $l
else 
	echo "$l" >> missing.txt
fi

done

wc -l missing.txt

#2.Find partially complete runs and delete such lines
awk '{print $3}' mmpbsa.txt | grep -E --line-number '^$' | cut -d: -f1 >> incomplete.txt

while read -r l
do
sed -n "$l"p mmpbsa.txt |  awk '{print $1}' >> missing.txt
done < incomplete.txt 
wc -l missing.txt
sed -i '/^$/d' missing.txt

wc -l missing.txt
END
#3. Find bad contacts
#cp uniq_mmpbsa.txt mmpbsa_duplicate.txt
#rm nmmpbsa.txt
column -t uniq_mmpbsa.txt > fmmpbsa.txt
sed -i "/^$/d" fmmpbsa.txt

grep .......................000 fmmpbsa.txt >> badcontacts.txt
sed -i "/.......................000/d" fmmpbsa.txt
sed -i "/nan/d" fmmpbsa.txt

#4. Correct bad contacts and sum the energies

#awk '{print $2}' badcontacts.txt | cut -b 1-10 >> fmmcorrected.txt

#paste badcontacts.txt fmmcorrected.txt | awk '{print $1 "\t" $3+$4+$5}' > fmmpbsa_cum.txt

awk '{print $1 "\t" $2+$3+$4+$5}' fmmpbsa.txt > mmpbsa_final.txt

