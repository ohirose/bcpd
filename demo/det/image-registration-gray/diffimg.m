disp('[3] Display results');

DIRBASE=sprintf('%s/../../..',     pwd);
DIRDATA=sprintf('%s/../../../data',pwd);

T=imread(sprintf('%s/111.jpg',DIRDATA));
S=imread(sprintf('%s/222.jpg',DIRDATA));
O=imread('L1.jpg');

imwrite(((S-T)+(T-S))*3,'diffdata.jpg');
imwrite(((O-T)+(T-O))*3,'diffL1.jpg');

O=imread('L2.jpg');
imwrite(((O-T)+(T-O))*3,'diffL2.jpg');

tile=imtile({'imy.jpg','imx.jpg','L1.jpg','L2.jpg','diffdata.jpg','diffL1.jpg','diffL2.jpg'},'GridSize',[1 7]);
imshow(tile);

