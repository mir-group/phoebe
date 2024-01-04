export QE_PATH=/path/to/qe/bin/

phono3py --qe -d --dim="3 3 3" -c scf.in
mkdir supercells
mv supercell* supercells
bash generate.sh

# HERE run the qe calculations!
# you can run them with more processes to speed this up
for i in $(seq -f "%05g" 1 49) ; do
  mpirun -np 4 $QE_PATH/pw.x -in disp-$i.in > disp-$i.out
done

# after running them, run these lines to collect the fc2 and fc3s
 phono3py --qe --cf3 disp-{00001..49}.out
phono3py --qe --dim="3 3 3" -c scf.in --sym-fc
