cd bcpd-nonrigid
sed -i.bak s/win=0\;/win=1\;/ *.m;
rm *.bak;
cd ..

cd bcpd-plusplus
sed -i.bak s/win=0\;/win=1\;/ *.m;
rm *.bak;
cd ..

cd bcpd-rigid
sed -i.bak s/win=0\;/win=1\;/ *.m;
rm *.bak;
cd ..
