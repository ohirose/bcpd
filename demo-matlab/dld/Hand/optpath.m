png=0;

fp=fopen('.optpath.bin');
N =fread(fp,  1,  'int32' );
D =fread(fp,  1,  'int32' );
M =fread(fp,  1,  'int32' );
L =fread(fp,  1,  'int32' );
T =fread(fp,D*M*L,'double');
X =fread(fp, D*N, 'double');
fclose(fp);

T =reshape(T,[D,M,L]);
X =reshape(X,[D,N]);

X=X';
bbox=[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))];

close
for l=1:L
  Y=T(:,:,l)';
  plot(X(:,1),X(:,2), 'bo','MarkerSize',8); hold on;
  plot(Y(:,1),Y(:,2),'-ro','MarkerSize',5,'MarkerFaceColor',[1,0,0]);
  axis(bbox);

  if(png==1)
    fn=sprintf('otw-%04d.png',l);
    print(fn,'-dpng');
  end;

  pause(0.1);
  hold off;
end

