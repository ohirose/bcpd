clc; clear; close all;
addpath(genpath('../../'));

DIRBASE=sprintf('%s/../../..',     pwd);
DIRDATA=sprintf('%s/../../../data',pwd);

fnm   =sprintf('%s/bcpd',        DIRBASE);
fnw   =sprintf('%s/win/bcpd.exe',DIRBASE);
if(ispc) EXE=fnw; else EXE=fnm; end;

ncol=30;

idy=1;
idx=3;

vx=sprintf('%s/tr_reg_%.3d.vert.txt',     DIRDATA,idx);
vy=sprintf('%s/tr_reg_%.3d.vert.txt',     DIRDATA,idy);
fx=sprintf('%s/tr_reg_%.3d.wks.c%.2d.txt',DIRDATA,idx,ncol);
fy=sprintf('%s/tr_reg_%.3d.wks.c%.2d.txt',DIRDATA,idy,ncol);
tr=sprintf('%s/tr_reg.face.txt',          DIRDATA);

IN =sprintf('-x%s -X%s -y%s -Y%s -Ggeo,.1,%s',vx,fx,vy,fy,tr);
OPT='-ux -Ux -n200 -c1e-6 -h -A -e.4 -l100 -b1 -w0.1 -g1 -r1 -sY -DB,2000,.05 ';
cmd=sprintf('%s %s %s',EXE,IN,OPT); system(cmd);

x=load(vx);
y=load(vy);
t=load('output_y.txt');

plot3(x(:,1),x(:,2),x(:,3),'.b'); hold on; 
plot3(t(:,1),t(:,2),t(:,3),'.r'); axis off; 
daspect([1 1 1]);

