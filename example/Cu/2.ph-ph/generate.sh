for i in $(seq -f "%05g" 1 49);
do
  cp template.in disp-$i.in
  cat supercells/supercell-$i.in >> disp-$i.in
done

