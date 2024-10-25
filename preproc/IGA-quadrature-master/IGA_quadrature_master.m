
function [xi,w]=IGA_quadrature_master(k,over,p,r)

nel=length(k)-1;
    n=(nel)*(p-r)+r+1;    
    [knot]=Sp_p_r(nel,n,p,r, k,over); 
    [w, xi, ~, ~] = getOptimalQuadPoints(knot, p );

end



function [knot]=Sp_p_r(nel,n,p,r,k,over)

knot(1:p+1)=k(1,1)*ones(1,p+1);
mid=ceil(nel/2);
ik=p+2;
if over
    %              knot(ik)=(k(1,2)+k(1,1))/2;
    %              ik=ik+1;
    knot(ik:ik+p-1)=k(1,1)+(k(1,2)-k(1,1))/(p+1)*(1:p);
    ik=ik+p;
end
if mid>1
    for i=1:mid-1
        for h=0:p-r-1
            knot(ik+h)=k(1,i+1);
        end
        ik=ik+p-r;
    end
else
    i = mid-1;
end
if mod(n,2)==1
    knot(ik)=(k(1,i+1)+k(1,i+2))/2;
    ik=ik+1;
end
for i=mid:nel-1
    for h=0:p-1-r
        knot(ik+h)=k(1,i+1);
    end
    ik=ik+p-r;
end
if over
    knot(ik:ik+p-1)=k(1,end-1)+(k(1,end)-k(1,end-1))/(p+1)*(1:p);
    ik=ik+p;
end
knot(ik:ik+p)=k(1,end)*ones(1,p+1);

end
