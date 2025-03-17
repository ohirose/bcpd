%disp('[3] Conversion: output text files => image');

function convOut(level)

  DIRBASE=sprintf('%s/../../..',pwd);
  DIRDATA=sprintf('%s/data',DIRBASE);
  
  fnx=sprintf('%s/111.jpg',DIRDATA);
  fny=sprintf('%s/222.jpg',DIRDATA);
  
  fx=imread(fnx);
  fy=imread(fny);
  
  [nc nr Df]=size(fx);
  Y=load(sprintf('L%d.txt',level));
  ft=zeros(nc,nr,Df);
  
  for i=1:nr
    for j=1:nc
      for d=1:Df
        idx=(j-1)+nc*(i-1)+1;
        ci=int32(Y(idx,1));
        cj=int32(Y(idx,2));
  
        if ci>=1 && cj>=1 && ci<=nr && cj<=nc 
          ft(ci,cj,d)=double(fy(i,j,d))/double(255);
        end
      end
    end
  end
  
  ft=medfilt2(ft,[5,5]);
  
  imwrite(fx,'imx.jpg');
  imwrite(fy,'imy.jpg');
  imwrite(ft,sprintf('L%d.jpg',level));
  
end

