cd /home/chemical/phd/chz198152/scratch/covid/pep-nsp7/simulations/7/
module load apps/gromacs/5.1.4/intel

name=md
gmx_mpi trjconv -f ${name}'.xtc' -s ${name}.tpr -o $name'cluster.gro' -e 0.001 -pbc cluster -n ind.ndx << EOF
1
0
EOF
gmx_mpi grompp -f ${name}.mdp -c ${name}cluster.gro -o ${name}cluster.tpr -p topol.top -n ind.ndx -maxwarn 2
gmx_mpi trjconv -f ${name}.xtc -o ${name}_cluster.xtc -s ${name}cluster.tpr -pbc nojump << EOF
0
EOF
gmx_mpi trjconv -f ${name}_cluster.gro -s ${name}cluster.tpr -n ind.ndx -o pept${name}_cluster.pdb  << EOF
1
EOF
END

rm -rf frames
mkdir frames

gmx_mpi trjconv -f md_cluster.xtc -o frames/frame.gro -s md.tpr -skip 50 -sep << EOF
1
EOF

for j in {1..200}
do
gmx_mpi sasa -f frames/frame$j.gro -s frames/frame$j.gro  -or frames/complex_r$j -oa frames/complex_a$j -o frames/complex_t$j -surface 1 -output 1 -ndots 2000 -n ind.ndx &

gmx_mpi sasa -f frames/frame$j.gro -s frames/frame$j.gro -or frames/nsp7_r$j -oa frames/nsp7_a$j -o frames/nsp7_t$j -surface 19 -output 19 -ndots 2000 -n ind.ndx &

gmx_mpi sasa -f frames/frame$j.gro -s frames/frame$j.gro -or frames/pep_r$j -oa frames/pep_a$j -o frames/pep_t$j -surface 20 -output 20 -ndots 2000 -n ind.ndx &

done
wait
cd ..
bash mmpbsa.sh
