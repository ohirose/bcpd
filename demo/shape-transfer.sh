Y=42
X=01

grep v ${Y}.obj|tr ' ' '\t'| cut -f1 > tmp1.txt; ## vertex symbols
grep f ${Y}.obj                      > tmp3.txt; ## face information
./bcpd -x shape01-rot180.txt -y shape42.txt                -J300 -K70 -p -r1 -c1e-6 -g 10 -n1000 -b2   -l2 -sY -DB,10000,0.08
./bcpd -x shape01-rot180.txt -y output_y.interpolated.txt  -J300 -K70 -p -r1 -c1e-6 -g.10 -n1000 -b1.2 -l2 -sY -uxd -o'rewind_'

for PSET in output_y.interpolated rewind_y; do
  paste tmp1.txt $PSET.txt |tr '\t' ' ' > tmp2.txt;
  cat tmp2.txt tmp3.txt > $PSET.obj
done;
rm tmp?.txt;
