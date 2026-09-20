close all; clear;

convIn;

disp('[2] Function registration'); pause(.5);
DIRBASE=sprintf('%s/../../..',     pwd);
DIRDATA=sprintf('%s/../../../data',pwd);

fnm   =sprintf('%s/bcpd',        DIRBASE);
fnw   =sprintf('%s/win/bcpd.exe',DIRBASE);
if(ispc) EXE=fnw; else EXE=fnm; end;

% Hierarchy 1
IN ='-x x.txt -y y.txt -X fx.txt -Y fy.txt';
OPT='-A -n500 -l20 -b.9 -w0.0 -c1e-6 -h -g3 -DB,10000,.01 -r1';
cmd=sprintf('%s %s %s',EXE,IN,OPT);
system(cmd);
copyfile output_y.txt L1.txt
convOut(1);

% Hierarchy 2
IN ='-x x.txt -y L1.txt -X fx.txt -Y fy.txt';
OPT='-A -n500 -l20  -b0.8 -w0.0 -c1e-6 -h -g.1 -DB,30000,.01 -r1';
cmd=sprintf('%s %s %s',EXE,IN,OPT);
system(cmd);
copyfile output_y.txt L2.txt
convOut(2);

diffimg;

