clear; close all;
%% input files
for i=1:40
  HAND_ID=i;
  X   =sprintf('%s/../../../data//hand%.3d.txt',pwd,HAND_ID);
  H   =sprintf('%s/Model/H-lout%.3d.txt',pwd,HAND_ID);
  mu  =sprintf('%s/Model/m-lout%.3d.txt',pwd,HAND_ID);
  
  DIRBASE=sprintf('%s/../../../',pwd);
  fnm =sprintf('%s/bcpd',        DIRBASE);
  fnw =sprintf('%s/win/bcpd.exe',DIRBASE);
  if(ispc); EXE=fnw; else; EXE=fnm; end
  
  %% parameters
  lmd ='1e-1';
  omg ='1e-2';
  gma ='1';
  
  %% execution
  cmd=sprintf('%s -un -x%s -y%s -C%s -l%s -w%s -g%s -h -sA -n50',EXE,X,mu,H,lmd,omg,gma);
  system(cmd); %optpath;

  %% parameters (BCPD)
  lmd ='10';
  bet ='.3';
  omg ='0';
  gma ='1';

  Y='output_y.txt';
  cmd=sprintf('%s -ux -x%s -y%s -l%s -b%s -w%s -g%s -h -n50 -sA',EXE,X,Y,lmd,bet,omg,gma);
  system(cmd); optpath;

end

