cd /home/chemical/phd/chz198152/scratch/covid/pep-nsp7/
module load apps/gromacs/5.1.4/intel

j=1
mkdir emposes

for i in {1..4000}
do

line=`sed "${i}q;d" names${j}.txt`

mkdir em$line
cd em$line
cp ../dockedposes/$line.pdb .
cp ../em.mdp .
cp ../ions.mdp .
cp -r ../charmm36-mar2019.ff .

gmx_mpi pdb2gmx -f $line.pdb -o dock$line.gro -p dock$line.top -ignh << EOF
1
1
EOF
gmx_mpi editconf -f dock$line.gro -o dock$line'_'n.gro -bt cubic -d 1.5 -c
gmx_mpi solvate -cp dock$line'_'n.gro -cs spc216.gro -p dock$line.top -o dock$line'_'s.gro
gmx_mpi grompp -f ions.mdp -c dock$line'_'s.gro -p dock$line.top -o ions$line.tpr 
gmx_mpi genion -s ions$line.tpr -o dock$line'_'i.gro -p dock$line.top -neutral -conc 0.15 << EOF
13
EOF
gmx_mpi grompp -f em.mdp -c dock$line'_'i.gro -p dock$line.top -o em$line.tpr
export OMP_NUM_THREADS=4
mpirun -np 6 gmx_mpi mdrun -v -s em$line.tpr -o em.trr -g em$line.log -c em$line.gro -ntomp 4

cp em$line.tpr em$line.log em$line.gro ../emposes/.

cd ..
rm -rf em$line 
done
