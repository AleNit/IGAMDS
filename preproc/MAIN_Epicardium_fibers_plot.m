
% visualize fiber orientation on the Epicardium surface

clc
clear 
close all
addpath('./IGA-quadrature-master/')


%% general input
np=5;                                               % number of patches
alp=pi/3;                                           % fiber orientation
col=0.8.*ones(1,3);


%% geometrical input
% surface
[p,q,U,V,CP]=readgeo1('./Epicardium_surface_NURBS.txt');


% plot initial geometry
figure(1)
el=false;                 % plot element boundaries
conp=false;               % plot control polygon
for ip=1:np
    plotNURBSsurf(p(ip),q(ip),U{ip},V{ip},CP{ip},el,conp,col)  
    create_el_edges2(p(ip),q(ip),U{ip},V{ip},CP{ip})
    hold on
end
% title('epicardium surface with fiber orientation','FontSize',14,'Interpreter','latex')


%% create fiber orientation
uu = linspace(0,1,11);
vv = uu;
for ip=1:np
    for i=1:length(uu)
        for j=1:length(vv)

            S0 = get_point_surf(p(ip),0,uu(i),U{ip},q(ip),0,vv(j),V{ip},CP{ip});

            [a1,a2]=epi_fibers(alp,p(ip),uu(i),U{ip},q(ip),vv(j),V{ip},CP{ip},ip);

            arrow3(S0,S0+a1'./1.8,'-r',0.6,1.5)     % fiber direction
            % arrow3(S0,S0+a2','-v',0.75,1.5)     % cross-fiber direction

        end
    end
end

axis([-4,4,-4,4,-8,4])
light('position',[60,60,60])
lighting gouraud
view(165,15)
