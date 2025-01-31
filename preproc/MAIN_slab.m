
% create input geomtry for the plane slab test

clc
clear 
close all
addpath('./IGA-quadrature-master/')


%% generate input
cname='test_planeslab';      % case name
np=1;                           % number of patches
wrt=true;                   % create input folder
cols={"#0072BD"};           % patch colors


%% geometrical input
% surface
p(1)=1; 
q(1)=1;
U{1}=[0 0 1 1];
V{1}=[0 0 1 1];
CP{1}(:,:,1)=[0.0 0.0; 2.0 2.0];
CP{1}(:,:,2)=[0.0,0.2; 0.0,0.2];
CP{1}(:,:,3)=[0.0,0.0; 0.0,0.0];
CP{1}(:,:,4)=[1.0,1.0; 1.0,1.0];

% patchwise integration parameters
over = 0; order = 4; regularity = 1;
disp(['in case of patchwise integration'])
disp(['estimated quadrature points per element: ',num2str(((order-regularity)/2)^2)])

% plot initial geometry
figure(1)
el=true;                 % plot element boundaries
conp=false;               % plot control polygon
for ip=1:np
    S1 = get_point_surf(p(ip),0,0,U{ip},q(ip),0,0,V{ip},CP{ip});
    S2 = get_point_surf(p(ip),0,0.01,U{ip},q(ip),0,0,V{ip},CP{ip});
    S3 = get_point_surf(p(ip),0,0,U{ip},q(ip),0,0.01,V{ip},CP{ip});
    plotNURBSsurf(p(ip),q(ip),U{ip},V{ip},CP{ip},el,conp,cols{ip})  
    hold on
    arrow3(S1,S1+(S2-S1)./norm(S2-S1)/3,'-r',1,2)
    arrow3(S1,S1+(S3-S1)./norm(S3-S1)/3,'-b',1,2)
end
title('initial surface','FontSize',14,'Interpreter','latex')



%% refinement
degu=2; 
degv=2;
refu=ceil(150/(length(U{1})-p(1)-2)); 
refv=ceil(3/(length(V{1})-q(1)-2));

ndof=0;
for ip=1:np

    [CP{ip},U{ip},V{ip},p(ip),q(ip)] = degree_elevate_surf(p(ip),q(ip),U{ip},V{ip},CP{ip},degu-p(ip),degv-q(ip));
    
    Ru = refinement_vec(U{ip},refu);
    Rv = refinement_vec(V{ip},refv);
    [CP{ip},U{ip},V{ip}] = knot_refine_surf(p(ip),q(ip),U{ip},V{ip},CP{ip},Ru,Rv);
    
    nu(ip) = length(CP{ip}(:,1,1));
    nv(ip) = length(CP{ip}(1,:,1));
    
    % compute quadrature points parameters
    [IPu{ip},IPv{ip},nqpu{ip},nqpv{ip}] = patchwise_int(U{ip},V{ip},CP{ip},over,order,regularity);

    ndofp(ip) = nu(ip)*nv(ip);
    ndof = ndof + ndofp(ip);

end

disp(['number of d.o.f.: ',num2str(ndof)])

% plot refined geometry
figure(2)
el=true;                 % plot element boundaries
conp=false;               % plot control polygon
for ip=1:np
    plotNURBSsurf(p(ip),q(ip),U{ip},V{ip},CP{ip},el,conp,cols{ip})  
end
title('refined surface','FontSize',14,'Interpreter','latex')



%% build connectivity matrix 
CONN = zeros(2,2);



%% write to output file
if (wrt)
    %---------------------------------------- create input folder
    infold=strcat('../',cname,'/input/');
    if (~exist(infold, 'dir'))
        mkdir(infold)
    end



    %---------------------------------------- model parameter: copy template
    LN={ ...
        '###################### CASE', ...
        strcat(num2str(np),'          		# number of patches (np)'), ...
        '.false.    		# read restart (is)', ...
        '200        		# step interval to print output (pint)', ...
        '100.0         	    # [ms] time interval for writing restart files (trest)', ...
        ' ', ...
        ' ', ...
        '###################### TISSUE PROPERTIES', ...
        '1			        # membrane model selector', ...
        '1.0        		# [mF/cm^3] tissue capacity (Cm)', ...
        '0.0001,0.0001  	# [S/cm] isotropic conductivity coefficient (Diso)', ...
        '1.0     		    # [1/cm] surface-to-volume cell ratio (chi)', ...
        ' ', ...
        ' ', ...
        '###################### DYNAMICS PARAMETERS', ...
        '0.008     		    # [ms] time-step size (dt)', ...
        '120.0        	    # [ms] runtime (rt)', ...
        '.false.			# reduced patchwise integration on reactive terms (redint)'};

    fid=fopen(strcat('../',cname,'/input/modelpar.in'),'w');
    for ii=1:length(LN)
        fprintf(fid,'%s \n',LN{ii});
    end
    fclose(fid);



    %---------------------------------------- stimulation protocol: copy template
    LM={ ...
        '############################ number of stimuli', ...
        '1', ...
        ' ', ...
        '############################ stimulus 1', ...
        '0.0                     	 	# t_start [ms]	', ...
        '0.5              		 	# t_end [ms]	', ...
        '1               			# patch ID', ...
        '0.0 0.0                       	 	# U subset', ...
        '0.0 0.0                   	 	# V subset', ...
        '0.0              			# Iapp [mA/cm^2]', ...
        '.true.                  		# overwrite potential', ...
        '1.0               			# value (dim.less or dimensional, depends on the membrane model)', ...
        '1,0,0,0                  		# forced side: {W,N,E,S}', ...
        ' ', ...
        '############################ stimulus 2' };
    
    % write to new file
    fid=fopen(strcat('../',cname,'/input/stim_prot.in'),'w');
    for ii=1:length(LM)
        fprintf(fid,'%s \n',LM{ii});
    end
    fclose(fid);



    %---------------------------------------- cell model: copy template / modify parameters
    filen=cellmodelname(str2num(LN{9}(1:2)));
    fid=fopen(filen);        
    FC=textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    
    % write to new file
    fid=fopen(strcat('../',cname,'/input/membrane.in'),'w');
    for ii=1:length(FC{1})
        fprintf(fid,'%s \n',FC{1}{ii});
    end
    fclose(fid);



    for ip=1:np

        %---------------------------------------- write geometry to file
        fname=strcat('../',cname,'/input/igeo_p',sprintf('%02d',ip),'.txt');    
        writesurf2file(fname,p(ip),q(ip),U{ip},V{ip},CP{ip})
    
    
    
        %---------------------------------------- quadrature point data to file
        fname=strcat('../',cname,'/input/qp_p',sprintf('%02d',ip),'.in');
        writeqp2file(fname,IPu{ip},IPv{ip},nqpu{ip},nqpv{ip})

    end


    %---------------------------------------- connectivity matrix for multipatch cases
    fileID=fopen(strcat('../',cname,'/input/conn_mp.in'),'w');    
    fprintf(fileID,'%i \n',length(CONN(:,1)));
    for i=1:length(CONN(:,1))
        fprintf(fileID,'%i \t %i \n',CONN(i,:));
    end
    fclose(fileID);


    %---------------------------------------- copy submission script
    dest=strcat('../',cname,'/go');
    status=copyfile('./template_files/go',dest);           

end





