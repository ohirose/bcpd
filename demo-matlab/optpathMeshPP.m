
% window size
w0=800*([0,0,1,1]);    % window (half): [left,bottom,width,height]
w1=800*([0,0,2,1]);    % window (full): [left,bottom,width,height]
% colors
ect='#0072BD';         % edge color (target)
ecs='k';               % edge color (source)
% figure titles
title1='Before Registration';
title2='Optimization Trajectory';
title3='After Registration';
t0=title1; w=w1; fg=figure('Name',t0,'NumberTitle','off'); set(fg,'Position',w);

% Rotate if FACE data
if tag==1
  vwc=[3,1,2];      % indicate view coordiates
  sgn=[-1,-1,1];    % indicate coordiate direction
  ang=(pi/180)*50;  % additional rotation around Y-axis
  % view-coordinates computation
  yc=[2,3,1];  % rotation around Y
  sgn=diag(sgn); 
  X0=X0(:,vwc)*sgn; 
  T0=T0(:,vwc)*sgn;
  T1=T1(:,vwc)*sgn;
  a=ang; Ry=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1]; Ry=Ry(yc,yc);
  if ang>0
    X0=X0*Ry';
    T0=T0*Ry';
    T1=T1*Ry';
  end
end

subplot(1,2,1); 
pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
pt=patch('vertices',T0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
axis equal; axis off;
title('Before Registration','FontSize',18);

subplot(1,2,2);
pt=patch('vertices',X0,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ect; hold on;
pt=patch('vertices',T1,'faces',f); pt.FaceAlpha=0; pt.EdgeColor=ecs; 
axis equal; axis off; 
title('After Registration','FontSize',18);

