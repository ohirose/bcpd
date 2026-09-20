close all; clear; clc;
X=load('Result/recovChef.txt');
K=load('Result/recovChefColor.txt');

N=size(X,1);
lim=300000;
mksz=8;
if N>lim; M=lim; else M=N; end
idx=randperm(N,M);
X=X(idx,:);
K=K(idx);


colors24 = [
    1.00, 0.00, 0.00;  % Red
    0.00, 1.00, 0.00;  % Green
    0.00, 0.00, 1.00;  % Blue
    1.00, 1.00, 0.00;  % Yellow
    1.00, 0.00, 1.00;  % Magenta
    0.00, 1.00, 1.00;  % Cyan
    0.50, 0.50, 0.50;  % Gray
    1.00, 0.50, 0.00;  % Orange
    0.50, 0.00, 0.50;  % Purple
    0.00, 0.50, 1.00;  % Sky Blue
    0.50, 1.00, 0.50;  % Light Green
    1.00, 0.75, 0.80;  % Pink
    0.63, 0.32, 0.18;  % Brown
    0.76, 0.69, 0.57;  % Tan
    0.93, 0.51, 0.93;  % Orchid
    0.50, 0.00, 0.00;  % Dark Red
    0.00, 0.50, 0.00;  % Dark Green
    0.00, 0.00, 0.50;  % Dark Blue
    0.50, 0.50, 0.00;  % Olive
    0.50, 0.00, 0.50;  % Dark Purple
    0.30, 0.30, 0.30;  % Dark Gray
    0.87, 0.72, 0.53;  % Light Brown
    0.69, 0.77, 0.87;  % Light Blue
    0.90, 0.91, 0.98   % Lavender
];

C=zeros(M,3);
for m=1:M
  C(m,:)=colors24(K(m),:);
end

scatter3(X(:,1),X(:,2),X(:,3),mksz,C,'filled');
daspect([1 1 1]);

