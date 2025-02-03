
function [a1,a2]=epi_fibers(alp,p,uu,U,q,vv,V,CP,ip)

g = base_vec(p,0,uu,U,q,0,vv,V,CP);    
if (ip==1 | ip==2 | ip==5) 
    a1 = g(:,1)/norm(g(:,1));
elseif (ip==3 | ip==4) 
    a1 = -g(:,1)/norm(g(:,1));
end
a2 = g(:,2)/norm(g(:,2));

% find normal array
n = cross(a1,a2);
n = n/sqrt(n'*n);

% rotate the first base vecor in the tangent plane until its z component becomes almost null
% use the Rodrigues' rotation formula
if ( a1(3) < 0 )
    th=-pi/500; 
else
    th=pi/500; 
end
while ( abs(a1(3)) > 1.0e-2 )
    a1 = a1.*cos(th) + cross(n,a1).*sin(th) + n.*(n'*a1)*(1-cos(th));                
end

% rotate again the base vector; rotation angle depend on the local surface orientation
th2=alp*abs(cos(n(3)));
a1 = a1.*cos(th2) + cross(n,a1).*sin(th2) + n.*(n'*a1)*(1-cos(th2));

% overwrite bitangent vector by a cross product
a2 = cross(a1,n);
a2 = a2/sqrt(a2'*a2);

end