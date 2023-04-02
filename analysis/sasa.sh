cd $PBS_O_WORKDIR
module load apps/gromacs/5.1.4/intel

rm -rf areas 
mkdir areas

for j in {1..17367}
do 

i=`sed "${j}q;d" x1.txt`

gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/complex_r${i} -o areas/complex_t${i} -surface 1 -output 1 -ndots 2000 -n ind1.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/nsp7_r${i} -o areas/nsp7_t${i} -surface 10 -output 10 -ndots 2000 -n ind1.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/pep_r${i} -o areas/pep_t${i} -surface 11 -output 11 -ndots 2000 -n ind1.ndx &


i=`sed "${j}q;d" x2.txt`

gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/complex_r${i} -o areas/complex_t${i} -surface 1 -output 1 -ndots 2000 -n ind2.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/nsp7_r${i} -o areas/nsp7_t${i} -surface 10 -output 10 -ndots 2000 -n ind2.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/pep_r${i} -o areas/pep_t${i} -surface 11 -output 11 -ndots 2000 -n ind2.ndx &
END
i=`sed "${j}q;d" x3.txt`

gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/complex_r${i} -o areas/complex_t${i} -surface 1 -output 1 -ndots 2000 -n ind3.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/nsp7_r${i} -o areas/nsp7_t${i} -surface 10 -output 10 -ndots 2000 -n ind3.ndx &
gmx_mpi sasa -f ../relax_dockedposes/outputs/${i}_0001.pdb -s ../relax_dockedposes/outputs/${i}_0001.pdb -or areas/pep_r${i} -o areas/pep_t${i} -surface 11 -output 11 -ndots 2000 -n ind3.ndx &

wait

rm areas/\#*

done

