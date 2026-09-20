clear; close all;
png=0; sub=1;

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

X=X'; Y0=T(:,:,1)';
bbox=[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))];

%% [left,bottom,width,height]
w0=800*([0,0,2,1]);
w1=800*([0,0,1,1]); pause(0.1);
if(sub==1) w1=w0; end
title1='Source and target point sets';
title2='Optimization trajectory';

%% figure 1
f1=figure('Name',title1,'NumberTitle','off'); set(f1,'Position',w1); figure(f1);
if(sub==1) subplot(1,2,1); end
plot(X (:,1),X (:,2),'bo','MarkerSize',8); hold on;
plot(Y0(:,1),Y0(:,2),'ro','MarkerSize',5,'MarkerFaceColor',[1,0,0]);
title(title1,'FontSize',18); axis off; hold off;

%% %% figure 2
if(sub==1) subplot(1,2,2); end
for l=1:L
  Y=T(:,:,l)';
  plot(X(:,1),X(:,2), 'bo','MarkerSize',8); hold on;
  plot(Y(:,1),Y(:,2), 'ro','MarkerSize',5,'MarkerFaceColor',[1,0,0]);
  title(title2,'FontSize',18); axis(bbox); axis off;

  if(png==1)
    fn=sprintf('Work/otw-%04d.png',l);
    print(fn,'-dpng');
  end

  pause(0.01);
  hold off;
end

