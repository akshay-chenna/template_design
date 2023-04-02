cd /home/chemical/phd/chz198152/scratch/covid/pep-nsp7/simulations/
module load apps/gromacs/5.1.4/intel
# Cat pbds generated by rosetta flexpepdock and find rmsd with the simulation


while read -r l
do
cd ../refined_new/$l/
rm all*
v=`echo $l | cut -d _ -f1`
cat *.pdb >> ../../simulations/allflexpep${v}.pdb

cd ..

done < ../refined_new/list.txt
