
% create input files for the Epicardium surface test

clc
clear 
close all
addpath('./IGA-quadrature-master/')


%% general input
cname='test_epi';                                   % case name
wrt=true;                                          % create input folder
np=5;                                               % number of patches
cols={"#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30"};


%% geometrical input
% surface
[p,q,U,V,CP]=readgeo1('./Epicardium_surface_NURBS.txt');

% patchwise integration parameters
over = 0; order = 4; regularity = 1;
disp(['in case of patchwise integration,'])
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
    arrow3(S1,S1+(S2-S1)./norm(S2-S1),'-r',1,2)
    arrow3(S1,S1+(S3-S1)./norm(S3-S1),'-b',1,2)
end
title('initial surface','FontSize',14,'Interpreter','latex')
lighting gouraud
drawnow



%% refinement
degu=2;     % same for all patches 
degv=2;
refu=ceil(180/(length(U{1})-p(1)-2)); 
refv=ceil(180/(length(V{1})-q(1)-2));

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


ndof=0;
for ip=1:np
    
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
drawnow



%% build connectivity matrix 
CONN=conn_mp_Epicardium(nu,nv,ndofp);



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
        '100        		# step interval to print output (pint)', ...
        '100.0         	    # [ms] time interval for writing restart files (trest)', ...
        ' ', ...
        ' ', ...
        '###################### TISSUE PROPERTIES', ...
        '3			        # membrane model selector', ...
        '1.4        		# [mF/cm^3] tissue capacity (Cm)', ...
        '0.0029,0.0014     		# [S/cm] isotropic conductivity coefficient (Diso)', ...
        '1400.0     		    # [1/cm] surface-to-volume cell ratio (chi)', ...
        ' ', ...
        ' ', ...
        '###################### DYNAMICS PARAMETERS', ...
        '0.025     		    # [ms] time-step size (dt)', ...
        '1000.0        	    # [ms] runtime (rt)', ...
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
        '2.5              		 	# t_end [ms]	', ...
        '1               			# patch ID', ...
        '0.47 0.53                       	 	# U subset', ...
        '0.37 0.43                   	 	# V subset', ...
        '0.025              			# Iapp [mA/cm^2]', ...
        '.false.                  		# overwrite potential', ...
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
