
function T=runbcpd(fnx,fny,prm,meth,fnf)
  bet=prm.bet;
  tau=prm.tau;
  gma=prm.gma;

  % binary
  fnm =sprintf('%s/../../bcpd',              pwd);
  fnw =sprintf('%s/../../win/bcpd.exe',      pwd);
  if(ispc) bcpd=fnw; else bcpd=fnm; end;

  %% common parameters
  omg ='0.0';
  lmd ='100';
  J   ='300';
  n   ='500';
  c   ='1e-6';
  nrm ='x';

  %% method parameters
  switch meth
    case 1 % bcpd
      tau='0';
      K  ='300';
      dwn='';
    case 2 % gbcpd
      K  ='300';
      dwn='';
    case 3 % gbcpd++
      K  ='300';
      dwn='-DB,3000,0.02';
  end
  prmm=sprintf('%s -Ggeo,%s,%s',dwn,tau,fnf);

  %% execution
  prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
  prm2=sprintf('-J%s -K%s -p -u%s',J,K,nrm);
  prm3=sprintf('-c%s -n%s -h -r1',c,n);
  cmd =sprintf('%s -x%s -y%s %s %s %s %s',bcpd,fnx,fny,prmm,prm1,prm2,prm3);
  system(cmd);
 
  %% output
  T=load('output_y.txt');
end 
