close all; clear;
addpath(genpath('../../'));

DIRBASE=sprintf('%s/../../..',pwd);
DIRDATA=sprintf('%s/data',DIRBASE);

fnm   =sprintf('%s/bcpd',        DIRBASE);
fnw   =sprintf('%s/win/bcpd.exe',DIRBASE);
if(ispc) EXE=fnw; else EXE=fnm; end;

vx=sprintf('%s/y.txt',  DIRDATA);
fx=sprintf('%s/fy2.txt',DIRDATA);
vy=sprintf('%s/skull.txt',DIRDATA);
fy=sprintf('%s/fy2.skull.txt',DIRDATA);

IN =sprintf('-x%s -X%s -y%s -Y%s',vx,fx,vy,fy);
OPT=sprintf('-ux -Tr -w0.75 -h -g3 -n500 -sY');

cmd=sprintf('%s %s %s',EXE,IN,OPT);
system(cmd); clear;

optpath3;
