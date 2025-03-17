close all; clear;
addpath('..');
%% input files
x   =sprintf('%s/../../data/chef_view001.txt',pwd);
y   =sprintf('%s/../../data/chef_view002.txt',pwd);
fnm =sprintf('%s/../../bcpd',                 pwd);
fnw =sprintf('%s/../../win/bcpd.exe',         pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.2';
gma ='1';
J   ='300';
f   ='0.3';
c   ='1e-6';
n   ='500';
dwn ='B,10000,0.04';
%% execution
prm1=sprintf('-w%s -g%s',omg,gma);
prm2=sprintf('-J%s -p -f%s -D%s',J,f,dwn);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -Tr -sT',bcpd,x,y,prm1,prm2,prm3);
system(cmd); 

Y=load(y);
X=load(x);
T=load('output_y.txt');

rigidresult
