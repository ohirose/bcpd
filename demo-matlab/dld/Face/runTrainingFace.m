clear; close all;
M=58; D=2; N=40; T=6;
K=50;

%% Create the 'Model' Directory
drc='Model';
if ~exist(drc,'dir')
  mkdir(drc);
end

X=zeros(D,M,T,N);
for n=1:N for t=1:T
  fnm=sprintf('../../../data/%.2d-%dm-opa.txt',n,t);
  fnf=sprintf('../../../data/%.2d-%df-opa.txt',n,t);
  if(exist(fnm)) fn=fnm; else fn=fnf; end;
  x=load(fn);
  X(:,:,t,n)=x';
end; end;
X=reshape(X,[D*M,N*T])'; 

% leave-one-person-out (leave six faces out)
for n=1:N
  del=(n-1)*T+1:n*T;
  Y=X; Y(del,:)=[]; %size(Y)
  L=size(Y,1);
  m=(sum(Y)/L);
  Y=Y-ones(L,1)*m;
  [B,Lmd]=eig(Y*Y');
  A=Y'*B;
  G=A*diag(1./sqrt(diag(A'*A)));
  Lmd=diag(Lmd/L);
  [dum,idx] = sort(Lmd,'descend');
  H=[Lmd G']'; H=H(:,idx); H=H(:,1:K);
  m=reshape(m,[D,M])';
  fnH=sprintf('Model/H-lout%.3d.txt',n); dlmwrite(fnH,H,'\t');
  fnm=sprintf('Model/m-lout%.3d.txt',n); dlmwrite(fnm,m,'\t');
end;

