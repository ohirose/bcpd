clear; close all;
%% input files
for i=1:40
  HAND_ID=i;
  DIRBASE=sprintf('%s/../../../',pwd);

  X   =sprintf('%s/data/hand%.3d.txt',   DIRBASE,HAND_ID);
  H   =sprintf('%s/Model/H-lout%.3d.txt',pwd,HAND_ID);
  mu  =sprintf('%s/Model/m-lout%.3d.txt',pwd,HAND_ID);
  
  fnm =sprintf('%s/bcpd',        DIRBASE);
  fnw =sprintf('%s/win/bcpd.exe',DIRBASE);
  if(ispc); EXE=fnw; else; EXE=fnm; end
  
  %% parameters
  omg ='1e-2';
  gma ='1';
  
  %% execution
  cmd=sprintf('%s -un -x%s -y%s -C%s -l1e-1 -w%s -g%s -h -sA -n50',EXE,X,mu,H,omg,gma);
  system(cmd);
  optpath;
end

