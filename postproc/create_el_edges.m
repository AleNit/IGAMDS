function create_el_edges(p,q,U,V,CP,str)
% draws the element edges
% J. Kiendl

mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

eps=10e-10;

% I) edges in u-direction
l=1;  % counting index of lines
for j2 = q+1:mv-q
  v = V(j2);
  j = findspan(v,V,nv);
  s = (U(mu-p)-U(p+1))/499;  %incremental step for u
  u = U(p+1);
  k=1;
  while u <= U(mu-p)+eps
    i = findspan(u,U,nu);
    Point(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    u=u+s;
  end
  l=l+1;
end

% II) edges in v-direction
for i2 = p+1:mu-p
  u = U(i2);
  i = findspan(u,U,nu);
  s = (V(mv-q)-V(q+1))/499;  %incremental step for v
  v = V(q+1);
  k=1;
  while v <= V(mv-q)+eps
    j = findspan(v,V,nv);
    Point(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP);
    k=k+1;
    v=v+s;
  end
  l=l+1;
end

if (str=='t')
    plot3 (Point(:,:,1),Point(:,:,2),Point(:,:,3),'--','Color','black');
elseif (str=='e')
    plot3 (Point(1:2:end,1:2:end,1),Point(1:2:end,1:2:end,2), ...
        Point(1:2:end,1:2:end,3),'-','Color','black');
else
    plot3 (Point(:,:,1),Point(:,:,2),Point(:,:,3),'Color','black');
end

end


