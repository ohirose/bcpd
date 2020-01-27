close all; clear;
%% set 'win=1' if windows
win=1;
%% input files
x   =sprintf('%s/../data/dragon-x.txt',pwd);
y   =sprintf('%s/../data/dragon-y.txt',pwd);
fnm =sprintf('%s/../bcpd',                pwd);
fnw =sprintf('%s/../win/bcpd.exe',        pwd);
if(win==1) bcpd=fnw; else bcpd=fnm; end;
%% parameters
omg ='0.0';
bet ='1.2';
lmd ='200';
gma ='10';
K   ='100';
J   ='300';
f   ='0.3';
c   ='1e-6';
n   ='500';
L   ='400';
nrm ='e';
dwn ='b,10000,0.02';
%% execution
prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
prm2=sprintf('-J%s -K%s -p -L%s -f%s -u%s -D%s',J,K,L,f,nrm,dwn);
prm3=sprintf('-c%s -n%s -h -r1',c,n);
cmd =sprintf('%s -x%s -y%s %s %s %s -sY',bcpd,x,y,prm1,prm2,prm3);
system(cmd);

X0=load(x);
T0=load(y);
T1=load('output_y.interpolated.txt');

% figure 2
f2=figure('Name','Before/After Registration','NumberTitle','off'); %set(f2,'Position',w1);
subplot(1,2,1);
plot3(T0(:,1),T0(:,2),T0(:,3),'.r','MarkerSize',3); hold on;
plot3(X0(:,1),X0(:,2),X0(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
title('Before Registration','FontSize',18);
subplot(1,2,2);
plot3(T1(:,1),T1(:,2),T1(:,3),'.r','MarkerSize',3); hold on;
plot3(X0(:,1),X0(:,2),X0(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
title('After Registration','FontSize',18);


