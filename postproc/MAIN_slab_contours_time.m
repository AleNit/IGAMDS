
% plot scalar contours over the simulated surface

clc
clear
close all
kf=1;



%% input parameters
loc='../test_planeslab/out/';         % case output folder
ns=0;                               % starting file id
nn=[200,200];                       % points on the surface for the contour representation
var=1;                              % variable to plot: 1-potential, 2-1st rec. variable, 3-iion, 4-iapp
wrf=false(1);                       % write frames to file



%% read some stuff from the input files
tmp=dir(strcat(loc,'p1_kin'));
ne=size(tmp,1)-4;

fname=strcat(loc,'../input/modelpar.in');
fileID=fopen(fname);
cac=textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);

splt=strsplit(cac{1}{2});
n_patches=str2num(splt{1});
splt=strsplit(cac{1}{4});
setpfig=str2num(splt{1});
splt=strsplit(cac{1}{16});
dt=str2num(splt{1});
dtfig=setpfig*dt;
clear cac


%% read geometry from input file
for np=1:n_patches

    filename=strcat(loc,'p',num2str(np),'_kin/refgeo.txt');    
    [p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np)]=read_refgeo(filename);
       
    [X(:,:,np),Y(:,:,np),Z(:,:,np)] = create_surf(p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np),nn);

end


figure(kf); kf=kf+1;
set(gcf,'position',[100,100,700,420])
for n=ns:1:ne

    t(n+1) = n*dtfig;    

    for np=1:n_patches

        cp_u=length(CP(:,1,1,np));
        cp_v=length(CP(1,:,1,np));
    
        % read solution file
        fname=strcat(loc,'p',num2str(np),'_kin/t',sprintf('%06d',n),'.h5');
        [CP_pot,CP_wrec,CP_ii,CP_ia]=getdatah5(fname,cp_u,cp_v);

        if (np==1)
            vv(n+1) = get_point_eval(p(np),0.3,U(:,np),q(np),0.5,V(:,np),CP(:,:,:,np),CP_pot);
        end
       
        % evaluate potential on a subgrid to get a smooth scalar field
        switch var
            case 1
                C = create_cont(p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np),CP_pot,nn);                
                titi='action potential [mV]';                
            case 2
                C = create_cont(p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np),CP_wrec(:,:,1),nn);
                titi='1st gating variable';
            case 3
                C = create_cont(p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np),CP_ii,nn);
                titi='ionic current [mA/cm^2]';
            case 4
                C = create_cont(p(np),q(np),U(:,np),V(:,np),CP(:,:,:,np),CP_ia,nn);
                titi='applied current [mA/cm^2]';
        end

        surf(X(:,:,np),Y(:,:,np),Z(:,:,np),C,'EdgeColor','none')
        hold on

    end

    colorbar
    % clim([-80,20])    
    axis equal
    % axis([0,2,0,0.2,-0.2,0.2])
    xlabel('x [cm]','FontSize',14,'Interpreter','latex');
    ylabel('y [cm]','FontSize',14,'Interpreter','latex');
    zlabel('z [cm]','FontSize',14,'Interpreter','latex');
    tit=strcat(titi,', t =',sprintf('%4.2f',n*dtfig),' [ms]'); 
    title(tit,'FontSize',14,'Interpreter','latex')
    drawnow   
    hold off    
    disp(['step ' num2str(n)])

    if (wrf)
        picname=strcat('./frames/frame_',sprintf('%05d',n));
        print('-dpng','-r200',picname);                        
    end    
    
end



%% compare data with https://doi.org/10.1016/j.cma.2016.12.022
vpeak=20;       % [mV]
vrest=-80;      % [mV]

fname='./planeslab_Goktepe_Khul.txt';
fileID = fopen(fname,'r');
cac = fscanf(fileID,'%f',[2 Inf])';
fclose(fileID);
tx=cac(:,1); vvx=cac(:,2);

figure(kf); kf=kf+1;
plot(tx,vvx.*(vpeak-vrest)+vrest,'o')
hold on
plot(t,vv,'-x')
xlabel('t [ms]','FontSize',14,'Interpreter','latex')
ylabel('v [mV]','FontSize',14,'Interpreter','latex')
legend('Gotkepe, Khull (2019)','present work')
