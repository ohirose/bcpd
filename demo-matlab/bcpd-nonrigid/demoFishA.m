close all; clear;
addpath('..');
%% input files
x   =sprintf('%s/../../data/fish-x.txt',pwd);
y   =sprintf('%s/../../data/fish-y.txt',pwd);
fnm =sprintf('%s/../../bcpd',           pwd);
fnw =sprintf('%s/../../win/bcpd.exe',   pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='2.0';
lmd ='2.0';
gma ='3.0';
c   ='1e-6';
n   ='500';
%% execution
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-c%s -n%s -h',c,n);
cmd =sprintf('%s -x%s -y%s %s %s -sA',bcpd,x,y,prm1,prm2);
system(cmd); optpath;
