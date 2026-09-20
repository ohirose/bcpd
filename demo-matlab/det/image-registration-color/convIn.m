close all; clear;

disp('[1] Conversion: input images => text files');

DIRBASE=sprintf('%s/../../..',pwd);
DIRDATA=sprintf('%s/data',DIRBASE);

fx=imread(sprintf('%s/antonio1-s.jpg',DIRDATA));
fy=imread(sprintf('%s/antonio2-s.jpg',DIRDATA));

[ncx nrx D]=size(fx);
[ncy nry D]=size(fy);
npx=ncx*nrx;
npy=ncy*nry;

cx=setprod([1:nrx],[1:ncx]);
cy=setprod([1:nry],[1:ncy]);

writematrix(cx,'x.txt','Delimiter','tab');
writematrix(cy,'y.txt','Delimiter','tab');

fx=reshape(fx,[npx D]);
fy=reshape(fy,[npy D]);

writematrix(fx,'fx.txt','Delimiter','tab');
writematrix(fy,'fy.txt','Delimiter','tab');
