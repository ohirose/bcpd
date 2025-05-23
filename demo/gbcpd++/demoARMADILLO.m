close all; clear; clc;
addpath('..');
%% input files
x   =sprintf('%s/../../data/armadillo-g-y.txt',pwd);
y   =sprintf('%s/../../data/armadillo-g-x.txt',pwd);
fnm =sprintf('%s/../../bcpd',        pwd);
fnw =sprintf('%s/../../win/bcpd.exe',pwd);
fnf =sprintf('%s/../../data/armadillo-g-triangles.txt',pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='1.0';
lmd ='50';
gma ='.1';
K   ='200';
J   ='300';
c   ='1e-6';
n   ='500';
L   ='100';
nrm ='x';
dwn ='B,10000,0.02';
tau ='1';
%% execution
kern=sprintf('geodesic,%s,%s',tau,fnf);
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-J%s -K%s -p -u%s -D%s',J,K,nrm,dwn);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -ux -G%s',bcpd,x,y,prm1,prm2,prm3,kern);
system(cmd);

X0=load(x);
T0=load(y);
T1=load('output_y.txt');
f =load(fnf); if min(min(f))==0; f=f+1; end;

tag=0;
optpathMeshPP;

