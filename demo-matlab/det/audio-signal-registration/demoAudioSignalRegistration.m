close all; clear;

DIRBASE=sprintf('%s/../../..',pwd);
DIRDATA=sprintf('%s/data',DIRBASE);

fnm   =sprintf('%s/bcpd',        DIRBASE);
fnw   =sprintf('%s/win/bcpd.exe',DIRBASE);
if(ispc) EXE=fnw; else EXE=fnm; end;

fn2=sprintf('%s/Ensoniq-ZR-76-Ac-Bass-2-C2.txt',          DIRDATA);
fn1=sprintf('%s/Alesis-Sanctuary-QCard-AcoustcBas-C2.txt',DIRDATA);

s2=load(fn2);
s1=load(fn1);

dlmwrite('t1.txt',[1:size(s1,1)]','precision','%d');
dlmwrite('t2.txt',[1:size(s2,1)]','precision','%d');
dlmwrite('s1.txt',s1,'precision','%d');
dlmwrite('s2.txt',s2,'precision','%d');

vx=sprintf('t2.txt'); t2=load(vx);
vy=sprintf('t1.txt'); t1=load(vy);
fx=sprintf('s2.txt'); s2=load(fx);
fy=sprintf('s1.txt'); s1=load(fy);

IN =sprintf('-x%s -X%s -y%s -Y%s',vx,fx,vy,fy);
OPT='-DB,2000,.05 -A -e.4 -l10 -b1 -w0.1 -g1 -r1 -ux -Ux -n1000 -c1e-6 -h';
cmd=sprintf('%s %s  %s',EXE,IN,OPT);
system(cmd);

o=load('output_y.txt');

t2=t2/44100;
t1=t1/44100;
o=o/44100;

f1=figure;
plot(t1,s1,'r'); hold on; plot(t2,s2,'b'); 
xlim([0,2.5]);
ylim([-1 1]);
set(gca,'FontSize',32); grid on; 
xlabel('Time (s)');
ylabel('Amplitude');

f2=figure;
plot(o,s1,'r'); hold on; plot(t2,s2,'b');
xlim([0,2.5]);
ylim([-1 1]);
set(gca,'FontSize',32); grid on;
xlabel('Time (s)');
ylabel('Amplitude');
