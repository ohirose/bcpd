disp('[3] Conversion: output text files => image');

DIRBASE=sprintf('%s/../../..',pwd);
DIRDATA=sprintf('%s/data',DIRBASE);

fnx=sprintf('%s/antonio1-s.jpg',DIRDATA);
fny=sprintf('%s/antonio2-s.jpg',DIRDATA);

fx=imread(fnx);
fy=imread(fny);

[nc nr Df]=size(fx);
Y=load('output_y.txt');
ft=zeros(nc,nr,Df);

for i=1:nr
  for j=1:nc
    for d=1:Df
      idx=(j-1)+nc*(i-1)+1;
      ci=int32(Y(idx,1));
      cj=int32(Y(idx,2));

      if ci>=1 && cj>=1 && ci<=nr && cj<=nc 
        ft(cj,ci,d)=double(fy(j,i,d))/double(255);
      end
    end
  end
end

fx=rot90(fx,3);
fy=rot90(fy,3);
ft=rot90(ft,3);

imwrite(fx,'imx.jpg');
imwrite(fy,'imy.jpg');
imwrite(ft,'imt.jpg');

tile=imtile({'imy.jpg','imx.jpg','imt.jpg'},'GridSize',[1 3]);
imshow(tile);

