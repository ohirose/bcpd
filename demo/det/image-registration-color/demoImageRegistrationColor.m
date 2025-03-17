close all; clear;

convIn;

disp('[2] Function registration'); pause(.5);
DIRBASE=sprintf('%s/../../..',pwd);

fnm   =sprintf('%s/bcpd',        DIRBASE);
fnw   =sprintf('%s/win/bcpd.exe',DIRBASE);
if(ispc) EXE=fnw; else EXE=fnm; end;

vx=sprintf('%s/x.txt', pwd);
vy=sprintf('%s/y.txt', pwd);
fx=sprintf('%s/fx.txt',pwd);
fy=sprintf('%s/fy.txt',pwd);

IN =sprintf('-x%s -X%s -y%s -Y%s',vx,fx,vy,fy);
OPT=sprintf('-A -l500 -w0.1 -b1.0 -h -g1 -n500 -DB,10000,.01 -r3');

cmd=sprintf('%s %s %s',EXE,IN,OPT);
system(cmd); clear;

convOut;
