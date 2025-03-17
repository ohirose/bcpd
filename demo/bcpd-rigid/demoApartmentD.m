close all; clear;
addpath('..');
%% input files
x   =sprintf('%s/../../data/apartment_004.txt',pwd);
y   =sprintf('%s/../../data/apartment_003.txt',pwd);
fnm =sprintf('%s/../../bcpd',                  pwd);
fnw =sprintf('%s/../../win/bcpd.exe',          pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.1';
gma ='0.1';
J   ='300';
e   ='0.3';
f   ='0.3';
c   ='1e-6';
n   ='500';
nrm ='n';
dwn ='B,2000,0.08';
%% execution
prm1=sprintf('-w%s -g%s',omg,gma);
prm2=sprintf('-J%s -p -e%s -f%s -u%s -D%s',J,e,f,nrm,dwn);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -Tr -sT',bcpd,x,y,prm1,prm2,prm3);
cmd
system(cmd); 

Y=load(y);
X=load(x);
T=load('output_y.txt');

rigidresult
