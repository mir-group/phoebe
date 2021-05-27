export QE_PATH="/your/path/to/QE/bin"

export NMPI=4
export NPOOL=4

ln -fs /your/path/to/thirdorder/thirdorder_espresso.py .

python3 thirdorder_espresso.py scf.in sow 2 2 2 -3 supercell_template.in

LIST_OF_FILES=DISP.supercell_template.in.*

for f in $LIST_OF_FILES; do
    mpirun -np $NMPI $QE_PATH/pw.x -npool $NPOOL -in $f > $f.out
done

find . -name 'DISP.supercell_template.in.*out' | sort -n | python3 thirdorder_espresso.py scf.in reap 2 2 2 -3
