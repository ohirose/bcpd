%clear; close all;
meshtype=0; png=0; sub=1; traj=1; rotb=0; rota=0; ttl=1;

%% read .optpath.bin
fp=fopen('.optpath.bin');
N =fread(fp,  1,  'int32' );
D =fread(fp,  1,  'int32' );
M =fread(fp,  1,  'int32' );
L =fread(fp,  1,  'int32' );
T =fread(fp,D*M*L,'double');
X =fread(fp, D*N, 'double');
T =reshape(T,[D,M,L]);
X =reshape(X,[D,N]);
% read mesh if available
[sz,ct]=fread(fp,2,'int32');
if ct==2 && sz(2)==3 
  nf=sz(1); nc=sz(2);
  f=fread(fp,nc*nf,'int32');
  f=reshape(f,[nc,nf])';
end
fclose(fp);
% node indices correction
if(min(min(f))==0) f=f+1; end

%% FACE
% view-coordinates computation
if tag==1
  vwc=[3,1,2];           % indicate view coordiates
  sgn=[-1,-1,1];         % indicate coordiate direction
  ang=(pi/180)*50;       % additional rotation around Y-axis
  yc=[2,3,1];  % rotation around Y
  sgn=diag(sgn); X=sgn*X(vwc,:); for l=1:L T(:,:,l)=sgn*T(vwc,:,l); end
  a=ang; Ry=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]; Ry=Ry(yc,yc);
  if ang>0
    X=Ry*X; for l=1:L T(:,:,l)=Ry*T(:,:,l); end
  end
end

%% view options
Lmax=200;              % view until l=Lmax
ndiv=100;              % #angles of the rotation view
ncam=1;                % #angles of the trajectory view (max=3);
w0=800*([0,0,1,1]);    % window (half): [left,bottom,width,height]
w1=800*([0,0,2,1]);    % window (full): [left,bottom,width,height]
ect='#0072BD';         % edge color (target)
ecs='k';               % edge color (source)
% figure titles
title1='Before Registration';
title2='Optimization Trajectory';
title3='After Registration';

%% bounding box computation
Y1=T(:,:,1)'; X1=X';
Yh=T(:,:,L)'; TT=reshape(T,[D,M*L])';
Xm=ones(N,  1)*mean(X1); Xc=X1-Xm;
Ym=ones(M,  1)*mean(X1); Yc=Y1-Ym;
Tm=ones(M*L,1)*mean(X1); Tc=TT-Tm;
Z=[Xc;Tc]; [V,E]=eig(Z'*Z); mxT=max(max(abs(Z*V)));
Z=[Xc;Yc]; [V,E]=eig(Z'*Z); mxY=max(max(abs(Z*V)));

%% rotation view: before
if rotb==1 
  t0=title1; w=w0; fg=figure('Name',t0,'NumberTitle','off'); set(fg,'Position',w);
  for i=0:ndiv
    a=2*pi*(i/ndiv); R=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1]; R=R(yc,yc);
    XR=X1*R; pt=patch('vertices',XR,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
    YR=Y1*R; pt=patch('vertices',YR,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
    if(ttl==1) title(t0,'FontSize',18); end
    axis(mxY*[-1 1 -1 1 -1 1]); axis off; axis vis3d; %axis equal;
    if png==1 
      fn=sprintf('Work/s0-%03d.png',i);
      %exportgraphics(fg,fn); % without margin
      print(fn,'-dpng');      % with margin
    end;
    pause(0.01); clf;
  end
  close all;
end 

%% optimization view
if traj==1 
  fg=figure('Name','Optimization Trajectory','NumberTitle','off');
  L =min(Lmax,L);
  nc=min(ncam,3); 
  for d=1:nc
    switch(d)
      case(1)
        R=eye(3); idx=[1 2 3];
      case(2)
        a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[3,1,2];
      case(3)
        a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[1,3,2];
    end
    Z=R(idx,idx)*X;
    for l=1:L U(:,:,l)=R(idx,idx)*T(:,:,l); end
    X0=Z'; Y0=U(:,:,1)';

    t0=sprintf('View %d',d);
    t1=sprintf('%s: View %d',title1,d);
    t2=sprintf('%s: View %d',title2,d);
    if sub==1  w=w1; else w=w0; end
    set(fg,'Position',w); figure(fg);

    for l=1:L
      Y=U(:,:,l)'; W=Z';

      if sub==1
        subplot(1,2,1);
        pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
        pt=patch('vertices',Y0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
        bbox=mxT*[-1 1 -1 1 -1 1]; axis equal; axis(bbox); axis off; 
        if ttl==1  title(t1,'FontSize',18); end
        subplot(1,2,2);
      end

      pt=patch('vertices',W,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
      pt=patch('vertices',Y,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
      bbox=mxT*[-1 1 -1 1 -1 1]; axis equal; axis(bbox); axis off; 
      if ttl==1  title(t2,'FontSize',18); end

      if png==1 
        axis square; axis tight; axis equal; axis off;
        fn=sprintf('Work/otw-%04d-v%d.png',l,d);
        %exportgraphics(fg,fn); % without margin
        print(fn,'-dpng');      % with margin
      end;

      pause(0.03); clf;
    end
  end
end
close all;

%% rotation view: after
if rota==1 
  t0=title3; w=w0; fg=figure('Name',t0,'NumberTitle','off'); set(fg,'Position',w);
  for i=0:ndiv
    a=2*pi*(i/ndiv); R=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1]; R=R(yc,yc);
    XR=X1*R; pt=patch('vertices',XR,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
    YR=Yh*R; pt=patch('vertices',YR,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
    axis(mxY*[-1 1 -1 1 -1 1]); axis off; axis equal;
    if ttl==1  title(t0,'FontSize',18); end
    if png==1 
      fn=sprintf('Work/s1-%03d.png',i);
      %exportgraphics(fg,fn); % without margin
      print(fn,'-dpng');      % with margin
    end;
    pause(0.03); clf;
  end
  close all;
end 

%% final view
if png==1 w=w0; else w=w1; end
fg=figure('Name','Before/After Registration','NumberTitle','off'); set(fg,'Position',w);
X0=X';

if meshtype==0
  if png==0 subplot(1,2,1); end; 
  pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
  pt=patch('vertices',Y1,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
  axis equal; axis off;
  if ttl==1 title('Before Registration','FontSize',18);  end
  if png==1 exportgraphics(fg,'shapes-before.png'); clf; end
  
  if png==0 subplot(1,2,2); end;
  pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
  pt=patch('vertices',Yh,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
  axis equal; axis off; 
  if ttl==1 title('After Registration','FontSize',18);  end
  if png==1 exportgraphics(fg,'shapes-after.png'); clf; end
else
  if png==0 subplot(1,2,1); end
  if meshtype==1
    C=ones(size(X0,1),1)*[1 1 1];
    pt=patch('vertices',X0,'faces',f,'FaceVertexCData',C);
    colormap gray(256); shading interp; lighting phong;
  else
    pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=1; pt.EdgeColor='k'; pt.FaceColor=[1 1 1];
  end
  camlight;
  axis equal; axis off;
  if ttl==1 title('Target Shape','FontSize',18); end
  if png==1 exportgraphics(fg,'shape-target.png'); clf; end
  
  if png==0 subplot(1,2,2); end
  if meshtype==1
    C=ones(size(X0,1),1)*[1 1 1];
    pt=patch('vertices',Yh,'faces',f,'FaceVertexCData',C);
    colormap gray(256); shading interp; lighting phong;
  else
    pt=patch('vertices',Yh,'faces',f); pt.FaceAlpha=1; pt.EdgeColor='k'; pt.FaceColor=[1 1 1];
  end
  camlight;
  axis equal; axis off; 
  if ttl==1 title('Deformed Shape','FontSize',18); end
  if png==1 exportgraphics(fg,'shape-deformed.png'); clf; end
end

 %% %% pcview
 %% if pcview==1 
 %%   Z=[Xc;Yc]; [V,E]=eig(Z'*Z);
 %%   [val,idx]=sort(diag(E),'descend'); V=V(:,idx);
 %%   X0=Xc*V'; 
 %%   Y0=Yc*V';  
 %%   for l=1:L
 %%    T(:,:,l)=V*T(:,:,l);
 %%   end
 %% end

