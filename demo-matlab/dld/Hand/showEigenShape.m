clear; close all; k=5; asp=-3;

m=load('Model/m-lout001.txt');
H=load('Model/H-lout001.txt');

Lmd=H(1,:)';
H=H(2:end,:);

h=sqrt(Lmd(k))*H(:,k); 
h=reshape(h,[56,2]);
s=m+asp*h; 

plot(s(:,1),s(:,2),'.-');
daspect([1 1 1])
