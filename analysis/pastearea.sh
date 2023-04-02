mapfile -t a < cent_struct.txt
<< 'END'
for i in {1..8}; do
sed -i '/^#/d' chC_r"${a[$i]}".xvg
sed -i '/^@/d' chC_r"${a[$i]}".xvg
awk '{print $2}' chC_r"${a[$i]}".xvg >> chC_r$i.txt 

sed -i '/^#/d' complex_r"${a[$i]}".xvg
sed -i '/^@/d' complex_r"${a[$i]}".xvg
tail -70 complex_r"${a[$i]}".xvg | awk '{print $2}' >> complexC_r$i.txt

tail -1 complex_t"${a[$i]}".xvg | awk '{ print $2 }' >> tot_$i.txt
done

paste chC_r{1..8}.txt >> monochC.txt
paste complexC_r{1..8}.txt >> complexC.txt
paste tot_{1..8}.txt >>tot.txt

END

for i in {1..8}; do
sed -i '/^#/d' chA_r"${a[$i]}".xvg
sed -i '/^@/d' chA_r"${a[$i]}".xvg
awk '{print $2}' chA_r"${a[$i]}".xvg >> chA_r$i.txt 

sed -i '/^#/d' complex_r"${a[$i]}".xvg
sed -i '/^@/d' complex_r"${a[$i]}".xvg
head -n -70 complex_r"${a[$i]}".xvg | awk '{print $2}' >> complexA_r$i.txt

done

paste chA_r{1..8}.txt >> monochA.txt
paste complexA_r{1..8}.txt >> complexA.txt


