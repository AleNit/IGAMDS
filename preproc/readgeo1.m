
% read multipatch geometry of the epicardium surface

function [p,q,U,V,CP]=readgeo1(fname)

np=5;       % number of patches
nc=6;       % number of control points per side

CP=cell(np,1);

for ip=1:np

    
    U{ip}=[0,0,0,0.25,0.5,0.75,1,1,1];
    V{ip}=U{ip};

    p(ip)=length(U{ip})-nc-1;
    q(ip)=p(ip);
    
    CP{ip}=zeros(nc,nc,4);
    
    fileID=fopen(fname);    
    cac=textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    
    if (isunix);    k=(ip-1)*nc^2+4;            % linux
    elseif (ispc);  k=(ip-1)*nc^2+2;    end     % windows

    for i=1:nc
        for j=1:nc
            tmp=strsplit(cac{1}{k});
            CP{ip}(i,j,1)=str2num(erase(tmp{1},["]",",","["]))*100;
            CP{ip}(i,j,2)=str2num(erase(tmp{2},["]",",","["]))*100;
            CP{ip}(i,j,3)=str2num(erase(tmp{3},["]",",","["]))*100;
            k=k+1;
        end
    end

    if (isunix);    k=(ip-1)*nc^2+186;            % linux
    elseif (ispc);  k=(ip-1)*nc^2+184;    end     % windows
    
    for i=1:6
        for j=1:6
            tmp=strsplit(cac{1}{k});
            CP{ip}(i,j,4)=str2num(tmp{1});
            k=k+1;
        end
    end

end


end