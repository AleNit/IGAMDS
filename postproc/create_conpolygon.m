function create_conpolygon(CP)
% plots control points and polygon
% J. Kiendl

% style
% lw - Line width
% ms - Marker size
% mf - Marker face coloring
% me - Marker edge coloring

% standard
% lw=0.5;
lw=1.5;
ms=6;
me='r';
mf='none';

% f?r dicken plot strip aus L-Platte
% lw=1.5;
% ms=6;
% me='k';
% mf='r';

% f?r dicke plots kragarm usw.
% lw=2.3;
% ms=8;
% me='k';
% mf='r';

% for k =1:length(CP(1,:,1))-1
%   for l =1:length(CP(:,1,1))-1
%     plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or','Linewidth',lw,...
%          'MarkerSize',ms,'MarkerEdgeColor',me,'MarkerFaceColor',mf);
%     plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or','Linewidth',lw,...
%          'MarkerSize',ms,'MarkerEdgeColor',me,'MarkerFaceColor',mf);
%   end
%   l=l+1;
%   plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or','Linewidth',lw,...
%          'MarkerSize',ms,'MarkerEdgeColor',me,'MarkerFaceColor',mf);
% end
% k=k+1;
% for l =1:length(CP(:,1,1))-1
%   plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or','Linewidth',lw,...
%          'MarkerSize',ms,'MarkerEdgeColor',me,'MarkerFaceColor',mf);
% end

% standard
for k =1:length(CP(1,:,1))-1
  for l =1:length(CP(:,1,1))-1
    plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or');
    plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or');
  end
  l=l+1;
  plot3(CP(l,k:k+1,1),CP(l,k:k+1,2),CP(l,k:k+1,3),'--or');
end
k=k+1;
for l =1:length(CP(:,1,1))-1
  plot3(CP(l:l+1,k,1),CP(l:l+1,k,2),CP(l:l+1,k,3),'--or');
end