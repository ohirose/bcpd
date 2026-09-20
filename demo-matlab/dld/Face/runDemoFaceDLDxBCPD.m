clear; close all; clc;

for HUMAN_ID=1:39
  for FACE_ID=1:6

    %% input files
    Xm  =sprintf('%s/../../../data/%.2d-%d%s-opa.txt',pwd,HUMAN_ID,FACE_ID,'m');
    Xf  =sprintf('%s/../../../data/%.2d-%d%s-opa.txt',pwd,HUMAN_ID,FACE_ID,'f');
    H   =sprintf('%s/Model/H-lout%.3d.txt',  pwd,HUMAN_ID);
    mu  =sprintf('%s/Model/m-lout%.3d.txt',  pwd,HUMAN_ID);
    
    if exist(Xm,'file')
      X=Xm; 
    else
      X=Xf; 
    end
    
    DIRBASE=sprintf('%s/../../../',pwd);
    fnm =sprintf('%s/bcpd',        DIRBASE);
    fnw =sprintf('%s/win/bcpd.exe',DIRBASE);
    if(ispc) EXE=fnw; else EXE=fnm; end
    
    %% parameters (DLD)
    lmd ='.5';
    omg ='1e-2';
    gma ='1';

    cmd=sprintf('%s -un -x%s -y%s -C%s -l%s -w%s -g%s -h -n50 -sA',EXE,X,mu,H,lmd,omg,gma)
    system(cmd); %optpath;

    %% parameters (BCPD)
    lmd ='10';
    bet ='.3';
    omg ='0';
    gma ='1';

    Y='output_y.txt'
    cmd=sprintf('%s -ux -x%s -y%s -l%s -b%s -w%s -g%s -h -n50 -sA',EXE,X,Y,lmd,bet,omg,gma)
    system(cmd); optpath;
  end
end
