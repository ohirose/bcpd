close all; clear;
addpath('..');
%% input files
x   =sprintf('%s/../../data/dragon-x.txt',pwd);
y   =sprintf('%s/../../data/dragon-y.txt',pwd);
fnm =sprintf('%s/../../bcpd',             pwd);
fnw =sprintf('%s/../../win/bcpd.exe',     pwd);
if(ispc) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='2.5';
lmd ='4.0';
gma ='10';
K   ='50';
J   ='300';
e   ='.25';
c   ='1e-6';
n   ='90';
%% execution
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-J%s -K%s -p',J,K);
prm3=sprintf('-c%s -n%s -e%s -h -r1 ',c,n,e);
cmd =sprintf('%s -x%s -y%s %s %s %s -sY',bcpd,x,y,prm1,prm2,prm3);
system(cmd); optpath3;

