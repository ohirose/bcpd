close all; clear;

disp('[1] Conversion: input images => text files');

DIRBASE=sprintf('%s/../../..',pwd);
DIRDATA=sprintf('%s/data',DIRBASE);

fx=imread(sprintf('%s/111.jpg',DIRDATA));
fy=imread(sprintf('%s/222.jpg',DIRDATA));

[nrx ncx]=size(fx);
[nry ncy]=size(fy);
npx=ncx*nrx;
npy=ncy*nry;

cx=setprod([1:nrx],[1:ncx]);
cy=setprod([1:nry],[1:ncy]);

writematrix(cx,'x.txt','Delimiter','tab');
writematrix(cy,'y.txt','Delimiter','tab');

fx=reshape(fx',[npx 1]);
fy=reshape(fy',[npy 1]);

writematrix(fx,'fx.txt','Delimiter','tab');
writematrix(fy,'fy.txt','Delimiter','tab');

