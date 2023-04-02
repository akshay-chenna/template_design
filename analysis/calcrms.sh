cut=0.4
while read -r line
do

for i in {1..5}
do

gmx rms -s $line/bbcapmodel$i.pdb -f ../templates/template$line.pdb -o $line'_'$i'.xvg' << EOF
0
0
EOF

a=`tail -1 $line'_'$i'.xvg' | awk '{ print $2 }'`

if (( $(bc <<<"$a <= $cut") ))
then
echo $line >> select.txt
fi

done

done < names.txt

uniq select.txt >> selected.txt
mkdir rmsd 
mv *.xvg rmsd/.
rm select.txt
cat selected.txt
wc -l selected.txt
