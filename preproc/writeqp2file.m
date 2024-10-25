
function writeqp2file(fname,IPu,IPv,nqpu,nqpv)

neu=length(nqpu);
nev=length(nqpv);
ngpu=sum(nqpu);
ngpv=sum(nqpv);



fID=fopen(fname,'w');

fprintf(fID,'%i \n',neu);
for i=1:neu
    fprintf(fID,'%i \n',nqpu(i));
end
fprintf(fID,'\n');
fprintf(fID,'%i \n',ngpu);
for i=1:ngpu
    fprintf(fID,'%16.14f \t %16.14f \n',IPu(i,1:2));
end
fprintf(fID,'\n');

fprintf(fID,'%i \n',nev);
for j=1:nev
    fprintf(fID,'%i \n',nqpv(j));
end
fprintf(fID,'\n');
fprintf(fID,'%i \n',ngpv);
for j=1:ngpv
    fprintf(fID,'%16.14f \t %16.14f \n',IPv(j,1:2));
end


fclose(fID);

end

