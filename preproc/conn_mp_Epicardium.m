
function CONNs=conn_mp_Epicardium(nu,nv,ndofp)

CONN = [];
c = 1;

%v edge [0:0,0:1] of patch 2 is slave of edge [1:1,0:1] of patch 1
sp=2;   % slave patch
mp=1;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*(i-1) + 1;
    master = nu(mp)*i;
    CONN(c,:) = [slave,master];
    c=c+1;
end

%v edge [1:1,0:1] of patch 3 is slave of edge [1:1,0:1] of patch 2
sp=3;   % slave patch
mp=2;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*i;
    master = sum(ndofp(1:mp-1)) + nu(mp)*i;
    CONN(c,:) = [slave,master];
    c=c+1;
end

%v edge [1:1,0:1] of patch 4 is slave of edge [0:0,0:1] of patch 3
sp=4;   % slave patch
mp=3;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*i;
    master = sum(ndofp(1:mp-1)) + nu(mp)*(i-1) + 1;
    CONN(c,:) = [slave,master];
    c=c+1;
end

%v edge [0:0,0:1] of patch 4 is slave of edge [0:0,0:1] of patch 1
sp=4;   % slave patch
mp=1;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*(i-1) + 1;
    master = nu(mp)*(i-1) + 1;
    CONN(c,:) = [slave,master];
    c=c+1;
end


%v edge [0:1,0:0] of patch 5 is slave of edge [0:1,1:1] of patch 1
sp=5;   % slave patch
mp=1;   % master patch
for i=1:nu(sp)
    slave = sum(ndofp(1:sp-1)) + i;
    master = (nv(mp)-1)*nu(mp) + i;
    CONN(c,:) = [slave,master];
    c=c+1;
end

%v edge [0:1,1:1] of patch 5 is slave of edge [0:1,1:1] of patch 3
sp=5;   % slave patch
mp=3;   % master patch
for i=1:nu(sp)
    slave = sum(ndofp(1:sp-1)) + (nv(sp)-1)*nu(sp) + i;
    master = sum(ndofp(1:mp-1)) + (nv(mp)-1)*nu(mp) + i;
    CONN(c,:) = [slave,master];
    c=c+1;
end

%v edge [0:0,0:1] of patch 5 is slave of edge [0:1,1:1] of patch 4
sp=5;   % slave patch
mp=4;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*(i-1) + 1;
    master = sum(ndofp(1:mp-1)) + (nv(mp)-1)*nu(mp) + i;
    CONN(c,:) = [slave,master];
    c=c+1;
end

% edge [1:1,0:1] of patch 5 is slave of edge [0:1,1:1] of patch 2
sp=5;   % slave patch
mp=2;   % master patch
for i=1:nv(sp)
    slave = sum(ndofp(1:sp-1)) + nu(sp)*i;
    master = sum(ndofp(1:mp-1)) + (nv(mp)-1)*nu(mp) + i;
    CONN(c,:) = [slave,master];
    c=c+1;
end


% special treatment for triple points; a slave node cannot be a master for
% another node! each slave must have only one master 
[vals,ic,jc] = unique(CONN(:,1),'stable');
valm = CONN(ic,2);

CONNs=[vals,valm];

id=find(CONNs(:,1)==sum(ndofp));
CONNs(id,2)=sum(ndofp(1:2));


end
