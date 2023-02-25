QE_PATH=/path/to/qe/bin

NMPI=4

for i in $(seq -f "%05g" 1 57) ; do
  echo $i
  mpirun -np $NMPI $QE_PATH/pw.x -in disp-$i.in > disp-$i.out
  # ideally, use regular QE and not phoebe-patched qe (if you've also done the elph tutorials)
  # if you however choose to, you need to remove the .out dir between runs to reset the gauge
  #rm -r out
done

