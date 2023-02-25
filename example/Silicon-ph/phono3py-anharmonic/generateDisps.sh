cp template.in disp.in
cat supercell.in >> disp.in
for i in $(seq -f "%05g" 1 57); do
  cp template.in disp-$i.in
  cat supercell-$i.in >> disp-$i.in
done
