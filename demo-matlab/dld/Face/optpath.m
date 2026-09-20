fp=fopen('.optpath.bin');
N =fread(fp,  1,  'int32' );
D =fread(fp,  1,  'int32' );
M =fread(fp,  1,  'int32' );
L =fread(fp,  1,  'int32' );
T =fread(fp,D*M*L,'double');
X =fread(fp, D*N, 'double');
fclose(fp);

T =reshape(T,[D,M,L]);
X =reshape(X,[D,N]);

X=X';
bbox=[min(X(:,1)),max(X(:,1)),min(X(:,2)),max(X(:,2))];

%% group definition for drawing
ng=7;
% mouse 
grp{1}=[ 1:13]; 
% eye
grp{2}=[14:21,14]; grp{3}=[22:29,22];
% eye brow
grp{4}=[30:34];    grp{5}=[35:39]; 
% mouse 
grp{6}=[40:47,40]; 
% nose
grp{7}=[48:58];


close
for l=1:L
  Y=T(:,:,l)';

  for g=1:ng
    plot(X(grp{g},1),X(grp{g},2),'-bo'); hold on;
    plot(Y(grp{g},1),Y(grp{g},2),'-ro','MarkerSize',5,'MarkerFaceColor',[1,0,0]);
    axis(bbox);
  end;

  pause(0.01);
  hold off;
end

