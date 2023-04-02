cd /home/chemical/phd/chz198152/scratch/covid/kb1-chC1/mmpbsa
module load apps/gromacs/5.1.4/intel

mmpbsa() {
		i=`sed "$1q;d" missing.txt`
		mkdir m${i}
		cd m${i}
		cp ../../emposes/em${i}.gro .
		cp ../residuetypes.dat .
		timeout 7m echo 11 10 | mpirun -np 1 -host `sed "$2q;d" $PBS_NODEFILE` .././g_mmpbsa -f em${i}.gro -s ../mmpbsa.tpr -n ../ind.ndx -i ../mmpbsa.mdp -pdie 4 -mme -mm mm.xvg -pbsa -pol pol.xvg -apol apol.xvg -decomp -mmcon contrib_mm_${i}.dat -pcon contrib_pol_${i}.dat -apcon contrib_apol_${i}.dat 


		vdw=`tail -n 1 mm.xvg | awk '{print $6}' `
		elec=`tail -n 1 mm.xvg | awk '{print $7}' `
		pol=`tail -n 1 pol.xvg | awk '{print $4-($3+$2)}'`
		apol=`tail -n 1 apol.xvg | awk '{print $4-($3+$2)}'`

		echo -e "${i}\t${vdw}\t${elec}\t${pol}\t${apol}" >> ../mmpbsa.txt
		
		cp contrib_mm_${i}.dat ../contrib/.
		cp contrib_pol_${i}.dat ../contrib/.
		cp contrib_apol_${i}.dat ../contrib/.
		
		cd ..
		rm -rf m${i}
}

for g in {1..1}; do
	for ((k=$(((g-1)*PBS_NTASKS+1)); k<=$(((g-1)*PBS_NTASKS+PBS_NTASKS)); k++)); do
		h=$((k%PBS_NTASKS+1))
		mmpbsa $k $h &
	done
	wait
done
