clear; close all;
png=0; sub=1; traj=1;

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

%% [left,bottom,width,height]
w0=800*([0,0,1,1]);
w1=800*([0,0,2,1]);
if(sub==1) w0=w1; end

title1='Before Registration';
title2='Optimization Trajectory';

if(traj==1)
  for d=1:D
    switch(d)
      case(1)
        R=eye(3); idx=[1,2,3];
      case(2)
        a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[3,1,2];
      case(3)
        a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[1,3,2];
    end
    Z=R(idx,idx)*X;
    for l=1:L U(:,:,l)=R(idx,idx)*T(:,:,l); end

    t0=sprintf('View %d',d);
    t1=sprintf('%s: View %d',title1,d);
    t2=sprintf('%s: View %d',title2,d);
    X0=Z'; Y0=U(:,:,1)';
    f1=figure('Name',t0,'NumberTitle','off'); set(f1,'Position',w0); figure(f1);

    if(sub==1)
      subplot(1,2,1);
      plot(X0(:,1),X0(:,2),'b.','MarkerSize',3); hold on;
      plot(Y0(:,1),Y0(:,2),'r.','MarkerSize',3,'MarkerFaceColor',[1,0,0]);
      daspect([1 1 1]); axis off; title(t1,'FontSize',18);
      subplot(1,2,2);
    end

    for l=1:L
      Y=U(:,:,l)'; W=Z';
  
      plot(W(:,1),W(:,2),'b.','MarkerSize',3); hold on;
      plot(Y(:,1),Y(:,2),'r.','MarkerSize',3,'MarkerFaceColor',[1,0,0]);
      bbox=[min(W(:,1)),max(W(:,1)),min(W(:,2)),max(W(:,2))];
      daspect([1 1 1]); axis(bbox); axis off; title(t2,'FontSize',18);

      if(png==1)
        fn=sprintf('Work/otw-%04d-v%d.png',l,d);
        print(fn,'-dpng');
      end;

      pause(0.01); hold off;
    end
  end
  close all;
end

X0=X';
T0=T(:,:,1)';
T1=T(:,:,L)';

% figure 2
f2=figure('Name','Before/After Registration','NumberTitle','off'); set(f2,'Position',w1);
subplot(1,2,1);
plot3(T0(:,1),T0(:,2),T0(:,3),'.r','MarkerSize',3); hold on;
plot3(X0(:,1),X0(:,2),X0(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
title('Before Registration','FontSize',18);
subplot(1,2,2);
plot3(T1(:,1),T1(:,2),T1(:,3),'.r','MarkerSize',3); hold on;
plot3(X0(:,1),X0(:,2),X0(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
title('After Registration','FontSize',18);

