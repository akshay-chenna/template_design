cd ~/scratch/covid/DESRES-Trajectory_sarscov2-10917618-no-water-zinc-glueCA/sarscov2-10917618-no-water-zinc-glueCA
module load apps/gromacs/2019.4/intel
gmx_mpi trjconv -f out.trr -s mae.pdb -o structures/o.pdb -sep << EOF
0
0
EOF

