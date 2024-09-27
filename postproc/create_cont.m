
function C=create_cont(p,q,U,V,CP,CPv,nn)

mu = length(U);
mv = length(V);
C=zeros(nn(1),nn(2));

% ploc=zeros(nn(1),nn(2),3);

eps=10e-10;

l=1;  % counting index of lines
r=(V(mv-q)-V(q+1))/(nn(2)-1);    %incremental step for v
v=V(q+1);
while v <= V(mv-q)+eps  
  s=(U(mu-p)-U(p+1))/(nn(1)-1);  %incremental step for u
  u=U(p+1);
  k=1;
  while u <= U(mu-p)+eps
    % ploc(k,l,:) = get_point_surf(p,0,u,U,q,0,v,V,CP);       
    vv = get_point_eval(p,u,U,q,v,V,CP,CPv);  
    C(k,l) = vv;
    k=k+1;
    u=u+s;
  end
  l=l+1;
  v=v+r;
end

end
