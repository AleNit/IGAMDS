
function [p,q,U,V,CP]=read_refgeo(filename)

fileID=fopen(filename);
formatspec='%s';
cac=textscan(fileID,formatspec,'Delimiter','\n');
fclose(fileID);

p=str2double(cac{1}(1));
q=str2double(cac{1}(2));
du=str2double(cac{1}(3));
dv=str2double(cac{1}(4));
cp_u=str2double(cac{1}(5));
cp_v=str2double(cac{1}(6));

U=zeros(du,1); V=zeros(dv,1);
CP=zeros(cp_u,cp_v,4);

for i=1:du
    U(i)=str2double(cac{1}(6+i));
end
for i=1:dv
    V(i)=str2double(cac{1}(6+du+i));
end
for i=1:4
    for j=1:cp_v
        for k=1:cp_u
            CP(k,j,i)=str2double(cac{1}(6+du+dv+k+(j-1)*cp_u+(i-1)*cp_v*cp_u));
        end
    end
end

end
