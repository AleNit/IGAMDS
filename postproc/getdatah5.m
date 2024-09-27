
% this assume a single current!

function [CP_pot,CP_wrec,CP_ii,CP_ia]=getdatah5(fname,cp_u,cp_v)

fid = H5F.open(fname);

dset_id = H5D.open(fid,'pot');
pot = H5D.read(dset_id);
H5D.close(dset_id);
 
dset_id = H5D.open(fid,'iion');
ii = H5D.read(dset_id);
H5D.close(dset_id);
 
dset_id = H5D.open(fid,'iapp');
ia = H5D.read(dset_id);
H5D.close(dset_id);

dset_id = H5D.open(fid,'wrec');
wrec = H5D.read(dset_id);
H5D.close(dset_id);
 
H5F.close(fid);


% assign fields to control point data structure
nw=length(wrec(1,:));
CP_pot=zeros(cp_u,cp_v);
CP_wrec=zeros(cp_u,cp_v,nw);
CP_ii=zeros(cp_u,cp_v);
CP_ia=zeros(cp_u,cp_v);

k=1;
for j = 1:cp_v
    for i = 1:cp_u
        CP_pot(i,j)=pot(k);        
        CP_wrec(i,j,:)=wrec(k,:);
        CP_ii(i,j)=ii(k);
        CP_ia(i,j)=ia(k);
        k=k+1;
    end
end

end
