#This is to check if any mmpbsa runs have been missed or need to be repeated

#cd $PBS_O_WORKDIR

#sort -u mmpbsa.txt > uniq_mmpbsa.txt
sed -i '/^$/d' bad_mmpbsa.txt

#1. Find completely  missed runs
e=`wc -l bad_mmpbsa.txt | awk '{print $1}'`
for ((j=1; j<=$((e));  j++))
do

l=`sed "${j}q;d" bad_mmpbsa.txt`

if grep -Fq "$l" bad_mmpbsa.txt
then
	echo $l >> bad_done.txt
else 
	echo "$l" >> bad_missing.txt
fi

done

wc -l bad_missing.txt

#2.Find partially complete runs and delete such lines
awk '{print $3}' bad_mmpbsa.txt | grep -E --line-number '^$' | cut -d: -f1 >> incomplete.txt

while read -r l
do
sed -n "$l"p bad_mmpbsa.txt |  awk '{print $1}' >> bad_missing.txt
done < incomplete.txt 
wc -l bad_missing.txt
sed -i '/^$/d' bad_missing.txt

wc -l bad_missing.txt

#3. Find bad contacts
grep .......................000 bad_mmpbsa.txt | awk '{ print $1}' > bad_contacts.txt
#column -t uniq_mmpbsa.txt > fmmpbsa.txt
#sed -i "/^$/d" fmmpbsa.txt

#grep .......................000 fmmpbsa.txt >> badcontacts.txt
#sed -i "/.......................000/d" fmmpbsa.txt
#sed -i "/nan/d" fmmpbsa.txt

#4. Correct bad contacts and sum the energies

#awk '{print $2}' badcontacts.txt | cut -b 1-10 >> fmmcorrected.txt

#paste badcontacts.txt fmmcorrected.txt | awk '{print $1 "\t" $3+$4+$5}' > fmmpbsa_cum.txt

#awk '{print $1 "\t" $2+$3+$4+$5}' fmmpbsa.txt > mmpbsa_final.txt

