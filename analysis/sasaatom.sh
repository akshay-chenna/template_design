cd /home/chemical/phd/chz198152/scratch/covid/nsp12-7-8/proteinonly
module load apps/gromacs/2019.4/intel

mkdir areas

for i in {0..8333}; 
do 

mpirun -np 1 gmx_mpi sasa -f structures/o$i.pdb -s structures/o$i.pdb -or areas/complex_r$i -o areas/complex_t$i -oa areas/complex_a$i -surface 20 -output 14 16 -ndots 2000 -n ind.ndx &
mpirun -np 1 gmx_mpi sasa -f structures/o$i.pdb -s structures/o$i.pdb -or areas/chA_r$i -o areas/chA_t$i -oa areas/chA_a$i -surface 14 -ndots 2000 -n ind.ndx &
mpirun -np 1 gmx_mpi sasa -f structures/o$i.pdb -s structures/o$i.pdb -or areas/chC_r$i -o areas/chC_t$i -oa areas/chC_a$i -surface 16 -ndots 2000 -n ind.ndx 
wait
done

