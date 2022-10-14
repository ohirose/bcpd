close all; clear;

a=2;
prm{1}.bet='1.2'; prm{1}.tau='0.5'; prm{1}.gma='3.0';
prm{2}.bet='0.7'; prm{2}.tau='0.0'; prm{2}.gma='0.1';
prm{3}.bet='0.3'; prm{3}.tau='0.0'; prm{3}.gma='0.1';

L=size(prm,2);

fny=sprintf('%s/../../data/female-y.txt',pwd);
fnx=sprintf('%s/../../data/female-x.txt',pwd);
fnf=sprintf('%s/../../data/female-triangles.txt',pwd);

% hierarchical registration
for l=1:L
  if l==1
    fnpsi=fny;
  else
    fnpsi=sprintf('%s/L%.2d.txt',pwd,l-1);
  end
  fnpso=sprintf('%s/L%.2d.txt',pwd,l);

  T=runbcpd(fnx,fnpsi,prm{l},a,fnf);
  dlmwrite(fnpso,T,'\t');
end

% visualization
X0=load(fnx);
T0=load(fny);
T1=load('output_y.txt');
f =load(fnf); if min(min(f))==0; f=f+1; end;

tag=0;
optpathMeshPP;

