type=0;

% figure 
f=figure('Name','Before/After Registration','NumberTitle','off');
if type==1
  subplot(1,2,1);
  plot3(Y(:,1),Y(:,2),Y(:,3),'.r','MarkerSize',3); hold on;
  plot3(X(:,1),X(:,2),X(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
  title('Before Registration','FontSize',18);
  subplot(1,2,2);
end

plot3(T(:,1),T(:,2),T(:,3),'.r','MarkerSize',3); hold on;
plot3(X(:,1),X(:,2),X(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
title('After Registration','FontSize',18);

