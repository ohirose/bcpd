Y=42;
X=01;
pfx1='transferA'
pfx2='transferB'
fn1=${pfx1}_y.interpolated;
fn2=${pfx2}_y;

grep v ${Y}.obj|tr ' ' '\t'| cut -f1  > v.txt;         ## vertex symbols
grep f ${Y}.obj                       > faces${Y}.txt; ## face information
grep v ${X}.obj|tr ' ' '\t'| cut -f2- > shape${X}.txt; ## vertices
grep v ${Y}.obj|tr ' ' '\t'| cut -f2- > shape${Y}.txt; ## faces
../../win//bcpd.exe -x shape${X}.txt -y shape${Y}.txt -h -J300 -K70 -p -r1 -L400 -c1e-6 -g 10 -n1000 -b2   -l2 -o${pfx1}_ -DB,10000,0.08
../../win//bcpd.exe -x shape${X}.txt -y ${fn1}.txt    -h -J300 -K70 -p -r1 -L400 -c1e-6 -g.10 -n1000 -b1.2 -l2 -o${pfx2}_ -ux 

for PSET in ${fn1} ${fn2}; do
  paste v.txt $PSET.txt |tr '\t' ' ' > tmp.txt;
  cat tmp.txt faces${Y}.txt > $PSET.obj
done;
rm tmp.txt shape??.txt v.txt faces??.txt

