close all; clear; clc;
X=load('Result/recovChef.txt');

N=size(X,1);
idx=randperm(N,300000);
X=X(idx,:);

colormap('turbo');
scatter3(X(:,1),X(:,2),X(:,3),4,X(:,3),'filled');
colorbar;
daspect([1 1 1]);

