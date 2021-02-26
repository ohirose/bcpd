close all; clear;
addpath('..');
%% set 'win=1' if windows
win=1;
%% input files
y   =sprintf('%s/../../data/face-x.txt',pwd);
x   =sprintf('%s/../../data/face-y.txt',pwd);
fnm =sprintf('%s/../../bcpd',           pwd);
fnw =sprintf('%s/../../win/bcpd.exe',   pwd);
if(win==1) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='0.3';
lmd ='1e4';
gma ='10';
K   ='150';
J   ='300';
c   ='1e-6';
n   ='500';
%% execution
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-J%s -K%s -p',J,K);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -sA',bcpd,x,y,prm1,prm2,prm3);
system(cmd); optpath3;

