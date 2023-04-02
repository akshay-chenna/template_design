while read -r line
do
cd $line

mkdir allmodels
tar -xf allModels.tar -C ./allmodels/
rm allModels.tar

for i in {1..5} 
do

sed -e '/TER/,$d' model$i.pdb >> modelpep$i.pdb
gmx pdb2gmx -f modelpep$i.pdb -o capmodel$i.pdb -ff charmm27 -water tip3p -ignh -ter << EOF
1
1
EOF
gmx trjconv -f capmodel$i.pdb -s capmodel$i.pdb -o bbcapmodel$i.pdb << EOF
4
EOF

done

rm *.top *.itp *.gro \#*

cd ..

done <  names.txt
