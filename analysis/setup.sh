cd /home/chemical/phd/chz198152/scratch/covid/pep-nsp7/simulations/5/
module load apps/gromacs/5.1.4/intel

name=top_5

gmx_mpi pdb2gmx -f ${name}.pdb -o conf.gro -p topol.top -ignh << EOF 
1
1
EOF
gmx_mpi editconf -f conf.gro -o complex_n.gro -bt cubic -d 1.5 -c
gmx_mpi solvate -cp complex_n.gro -cs spc216.gro -p topol.top -o complex_s.gro
gmx_mpi grompp -f ions.mdp -c complex_s.gro -p topol.top -o ions.tpr
gmx_mpi genion -s ions.tpr -o complex_i.gro -p topol.top -neutral -conc 0.15 << EOF
13
EOF
echo 'a 1-1116' >> createind.txt
echo a 1117-` wc -l conf.gro | awk '{print $1-3}'` >> createind.txt
echo name 19 NSP7 >> createind.txt
echo name 20 PEP >> createind.txt
echo '19 & 3' >> createind.txt
echo '20 & 3' >> createind.txt
echo q >> createind.txt
gmx_mpi make_ndx -f complex_i.gro -o ind.ndx < createind.txt
gmx_mpi genrestr -f complex_i.gro -n ind.ndx -o cansp7 << EOF 
21
EOF
cat addcaposrensp7.txt >> topol_Protein_chain_A.itp
echo '[ POSRE_CAPEP ]' >> ind.ndx
tail -n +1119 conf.gro | grep -n CA | cut -d : -f1 | tr '[\n]' ['\t'] >> ind.ndx
gmx_mpi genrestr -f complex_i.gro -n ind.ndx -o capep << EOF
23
EOF
cat addcaposrepep.txt >> topol_Protein_chain_B.itp

rm \#*




