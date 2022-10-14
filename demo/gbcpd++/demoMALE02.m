close all; clear; clc;
addpath('..');
%% input files
x   =sprintf('%s/../../data/male-y.txt',   pwd);
y   =sprintf('%s/../../data/male-x.txt',   pwd);
fnm =sprintf('%s/../../bcpd',              pwd);
fnw =sprintf('%s/../../win/bcpd.exe',      pwd);
fnf =sprintf('%s/../../data/male-triangles.txt',pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='1.3';
lmd ='100';
gma ='3';
K   ='300';
J   ='300';
c   ='1e-6';
n   ='500';
nrm ='x';
dwn ='B,3000,0.02';
tau ='.5';
%% execution
kern=sprintf('geodesic,%s,%s',tau,fnf);
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-J%s -K%s -p -u%s -D%s',J,K,nrm,dwn);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -ux -G%s',bcpd,x,y,prm1,prm2,prm3,kern);
system(cmd);

X0=load(x);
T0=load(y);
T1=load('output_y.interpolated.txt');
f =load(fnf); if min(min(f))==0; f=f+1; end;

tag=0;
optpathMeshPP;

