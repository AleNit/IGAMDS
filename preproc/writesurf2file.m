
function writesurf2file(fname,p,q,U,V,CP)

cpu=length(CP(:,1,1));
cpv=length(CP(1,:,1));

fID=fopen(fname,'w');

fprintf(fID,'#-------------------- U-DIR INPUT : p, du, [U] \n');
fprintf(fID,'%i \n',p);
fprintf(fID,'%i \n',length(U));
fprintf(fID,'%16.14f ',U);
fprintf(fID,'\n \n');

fprintf(fID,'#-------------------- V-DIR INPUT : q, dv, [V] \n');
fprintf(fID,'%i \n',q);
fprintf(fID,'%i \n',length(V));
fprintf(fID,'%16.14f ',V);
fprintf(fID,'\n \n');

fprintf(fID,'%i %i \n',[cpu,cpv]);
for i=1:cpu
    for j=1:cpv
        fprintf(fID,'%16.14f ',CP(i,j,:));
        fprintf(fID,'\n');
    end
end



fclose(fID);

end