cd /home/chemical/phd/chz198152/scratch/covid/nsp12-7-8/proteinonly/areas


for i in {0..8333}; do

sed -i '/^#/d' complex_r$i.xvg
sed -i '/^#/d' chA_r$i.xvg
sed -i '/^#/d' chC_r$i.xvg

sed -i '/^@/d' complex_r$i.xvg
sed -i '/^@/d' chA_r$i.xvg
sed -i '/^@/d' chC_r$i.xvg

awk '{ print $2 }' complex_r$i.xvg >> bound_r$i.txt
awk '{ print $2 }' chA_r$i.xvg >> unbound_r$i.txt
awk '{ print $2 }' chC_r$i.xvg >> unbound_r$i.txt
paste bound_r$i.txt unbound_r$i.txt >> temp.txt
awk '{print $1-$2}' temp.txt >> delsasa_r$i.txt
rm temp.txt

done

#paste delsasa_r{0..8333}.txt >> delsasa.txt
