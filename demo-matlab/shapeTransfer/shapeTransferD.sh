Y=09;
X=01;
pfx1='transferV1'
pfx2='transferV2'
fn1=${pfx1}_y;
fn2=${pfx2}_y;

STR="$(uname -s)"
case "${STR}" in
    MINGW*) ENV=MINGW;;
    *)      ENV=OTHERS;;
esac

if [[ ${ENV} == MINGW ]]
then
  EXE='../../win/bcpd.exe -h';
else
  EXE=../../bcpd;
fi;

grep v ${Y}.obj|tr ' ' '\t'| cut -f1  > v.txt;         ## vertex symbols
grep f ${Y}.obj                       > faces${Y}.txt; ## face information
grep v ${X}.obj|tr ' ' '\t'| cut -f2- > shape${X}.txt; ## vertices
grep v ${Y}.obj|tr ' ' '\t'| cut -f2- > shape${Y}.txt; ## faces
$EXE -x shape${X}.txt -y shape${Y}.txt -J300 -K70 -p -r1 -c1e-6 -g 10 -n1000 -b2   -l50 -o${pfx1}_ -DB,4000,0.08
$EXE -x shape${X}.txt -y ${fn1}.txt    -J300 -K70 -p -r1 -c1e-6 -g.10 -n1000 -b1.2 -l50 -o${pfx2}_ -DB,4000,0.08 -ux

for PSET in ${fn1} ${fn2}; do
  paste v.txt $PSET.txt |tr '\t' ' ' > tmp.txt;
  cat tmp.txt faces${Y}.txt > $PSET.obj
done;
rm tmp.txt shape??.txt v.txt faces??.txt transfer*.txt

