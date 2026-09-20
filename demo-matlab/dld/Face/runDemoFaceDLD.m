clear; close all;

%% parameters
omg ='1e-2';
gma ='1';

for HUMAN_ID=1:40
  for FACE_ID=1:6

    %% input files
    Xm  =sprintf('%s/../../../data/%.2d-%d%s-opa.txt',pwd,HUMAN_ID,FACE_ID,'m');
    Xf  =sprintf('%s/../../../data/%.2d-%d%s-opa.txt',pwd,HUMAN_ID,FACE_ID,'f');
    H   =sprintf('%s/Model/H-lout%.3d.txt',  pwd,HUMAN_ID);
    mu  =sprintf('%s/Model/m-lout%.3d.txt',  pwd,HUMAN_ID);
    
    if exist(Xm,'file'); X=Xm; else X=Xf; end
    
    DIRBASE=sprintf('%s/../../../',pwd);
    fnm =sprintf('%s/bcpd',        DIRBASE);
    fnw =sprintf('%s/win/bcpd.exe',DIRBASE);
    if(ispc) EXE=fnw; else EXE=fnm; end
    
    %% execution
    cmd=sprintf('%s -un -x%s -y%s -C%s -l1 -w%s -g%s -h -sA -n50',EXE,X,mu,H,omg,gma)
    system(cmd); 
    optpath;
  end
end
