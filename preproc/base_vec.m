function g = base_vec(p,i,u,U,q,j,v,V,CP)
% returns the base vectors g=[g1,g2] at point u,v
% i,j are the knot span indices. if unknown, insert 0!
% J. Kiendl

if (i==0); i = findspan(u,U,length(CP(:,1,1)));  end
if (j==0); j = findspan(v,V,length(CP(1,:,1)));  end

%  derive basis functions dR/du and dR/dv
N = deriv(i,p,u,U);    % basisfunc in u
M = deriv(j,q,v,V);    % basisfunc in v

sum = 0;
dsum = zeros(2,1);

for c = 0:q
  for b = 0:p
    R(b+1,c+1) = N(1,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    sum = sum + R(b+1,c+1);
            
    % derivatives
    dRl(b+1,c+1,1) = N(2,b+1)*M(1,c+1)*CP(i-p+b,j-q+c,4);
    dsum(1)  = dsum(1) + dRl(b+1,c+1,1);
    dRl(b+1,c+1,2) = N(1,b+1)*M(2,c+1)*CP(i-p+b,j-q+c,4);
    dsum(2)  = dsum(2) + dRl(b+1,c+1,2);
  end
end

    % divide through by sum
for c = 0:q
  for b = 0:p
    dRl(b+1,c+1,1) = dRl(b+1,c+1,1)/sum - (R(b+1,c+1)*dsum(1))/(sum^2);
    dRl(b+1,c+1,2) = dRl(b+1,c+1,2)/sum - (R(b+1,c+1)*dsum(2))/(sum^2);
  end
end

g = zeros(3,2);
for c = 0:q 
  for b = 0:p
    for a = 1:3
      g(a,1) = dRl(b+1,c+1,1)*CP(i-p+b,j-q+c,a) + g(a,1);
      g(a,2) = dRl(b+1,c+1,2)*CP(i-p+b,j-q+c,a) + g(a,2);
    end
  end
end

end