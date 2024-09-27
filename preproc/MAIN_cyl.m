
% create input geomtry for the cylindrical surface test

clc
clear 
close all


%% input parameters
cname='test_cyl';      % case name
wrt=true(1);                   % create input folder
fname=strcat('../',cname,'/input/igeo_p1.txt');   % gemetry file name

% initial geometry
R=1/sqrt(2);        % base radius
th=1/2*pi;    % base angle
L=1;        % cylinder height

p=2; 
q=1;
U=[0 0 0 1 1 1];
V=[0 0 1 1 ];
CP(:,:,1)=[0, 0; R*sin(th/2), R*sin(th/2); 2*R*sin(th/2), 2*R*sin(th/2)];
CP(:,:,2)=[0, L; 0, L; 0, L];
CP(:,:,3)=[0.0, 0.0; R*sin(th/2)*tan(th/2), R*sin(th/2)*tan(th/2); 0.0, 0.0];
CP(:,:,4)=[1.0, 1.0; cos(th/2), cos(th/2); 1.0, 1.0];

figure(1)
el=true(1);                 % plot element boundaries
conp=true(1);               % plot control polygon
plotNURBSsurf(p,q,U,V,CP,el,conp)
title('initial surface','FontSize',14,'Interpreter','latex')



%% refinement
degu=2; degv=2;
refu=120; refv=120;

[CP,U,V,p,q] = degree_elevate_surf(p,q,U,V,CP,degu-p,degv-q);

Ru = refinement_vec(U,refu);
Rv = refinement_vec(V,refv);
[CP,U,V] = knot_refine_surf(p,q,U,V,CP,Ru,Rv);

nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

figure(2)
el=true(1);                 % plot element boundaries
conp=false(1);               % plot control polygon
plotNURBSsurf(p,q,U,V,CP,el,conp)
title('refined surface','FontSize',14,'Interpreter','latex')


%% write to output file
if (wrt)

    % create input folder
    infold=strcat('../',cname,'/input/');
    if (~exist(infold, 'dir'))
        mkdir(infold)
    end

    % write geometry to file
    writesurf2file(fname,p,q,U,V,CP)

    % copy template input files
    dest=strcat('../',cname,'/go');
    status=copyfile('./template_files/go',dest);           

    dest=strcat('../',cname,'/input/modelpar.in');
    status=copyfile('./template_files/modelpar.in',dest); 

    dest=strcat('../',cname,'/input/stim_prot.in');
    status=copyfile('./template_files/stim_prot.in',dest); 

end


