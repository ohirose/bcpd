%% The MATLAB code trains hand shape models in a leave-one-out cross-validation manner.
clear; close all;

% drc: Directory to save results ('Model').
% L:   Total number of input files.
% K:   Number of eigenshapes to retain.

drc='Model';L=40; K=38;

%% Create the 'Model' Directory
if ~exist(drc, 'dir')
    mkdir(drc);
end

%% Reading Input Files 
% The following loop iterates through L files (hand001.txt to hand040.txt in the Data directory).
% Each file's data is loaded into the Xorig array, which stores shape data. The assumption is 
% that the points are the files are in a consistent format.
 
ct=0;
for l=1:L
  x=load(sprintf('../../../data/hand%.3d.txt',l));
  if(ct==0) [M,D]=size(x); Xorig=zeros(D,M,L); ct=ct+1; end;
  Xorig(:,:,l)=x'; 
end;

%% Shape Model Construction (for leave-one-out cross-validation) 
% A shape model is constructed for each dataset, from which the current shape data 
% (corresponding to l) is excluded. 
% m:   Mean shape   (M x D) 
% H:   Shape variations, defined as H=[Lmd Q']' 
% Lmd: Eigenvalues  (K)
% Q:   Eigenvectors (DM x K)

for l=1:L
  sprintf('l=%.3d\n',l)
  %% leave-one-out
  X=Xorig; X(:,:,l)=[];
  
  %% computing eigenshapes
  X=reshape(X,[D*M,L-1])'; 
  m=(sum(X)/(L-1));
  X=X-ones(L-1,1)*m;
  [B,Lmd]=eig(X*X');
  A=X'*B;
  Q=A*diag(1./sqrt(diag(A'*A)));
  Lmd=diag(Lmd/(L-1)+1e-10);
  [dum,idx] = sort(Lmd,'descend');
  H=[Lmd Q']'; H=H(:,idx); H=H(:,1:K);
  m=reshape(m,[D,M])';
  
  writematrix(H,sprintf('%s/H-lout%.3d.txt',drc,l),'Delimiter','tab');
  writematrix(m,sprintf('%s/m-lout%.3d.txt',drc,l),'Delimiter','tab');
end;

