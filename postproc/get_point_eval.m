
function vv = get_point_eval(p,u,U,q,v,V,CP,CPv)

i = findspan(u,U,length(CP(:,1,1)));
j = findspan(v,V,length(CP(1,:,1)));

Nu=basisfunc(i,p,u,U);
Nv=basisfunc(j,q,v,V);

SumNw = 0;

for c = 0:q
  for b = 0:p
    SumNw = Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)+SumNw;
  end
end

k=0;
R=zeros((q+1)*(p+1),1);
for c=0:q
    for b=0:p
        k=k+1;
        R(k)=Nu(b+1)*Nv(c+1)*CP(i-p+b,j-q+c,4)/SumNw;
    end
end

vv=0.0;
k=0;
for c=0:q
    for b=0:p
        k=k+1;
        vv=vv+R(k)*CPv(i-p+b,j-q+c);
    end
end

end