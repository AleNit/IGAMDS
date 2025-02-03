
% plot scalar contours over the simulated surface

clc
clear
close all
kf=1;



%% input parameters
loc='../test_epi/out/';         % case output folder
ns=0;                               % starting file id
nn=[200,200];                       % points on the surface for the contour representation
var=1;                              % variable to plot: 1-potential, 2-1st rec. variable, 3-iion, 4-iapp
wrf=false;                       % write frames to file



%% read some stuff from the input files
tmp=dir(strcat(loc,'p01_kin'));
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
for ip=1:n_patches

    filename=strcat(loc,'p',sprintf('%02d',ip),'_kin/refgeo.txt');    
    [p(ip),q(ip),U{ip},V{ip},CP{ip}]=read_refgeo(filename);
       
    [X(:,:,ip),Y(:,:,ip),Z(:,:,ip)] = create_surf(p(ip),q(ip),U{ip},V{ip},CP{ip},nn);

end


figure(kf); kf=kf+1;
set(gcf,'position',[100,100,1000,450])

TT=tiledlayout(1,2);
vv(1)=-83;

for n=1:1:ne

    t(n+1) = n*dtfig;    

    for ip=1:n_patches

        cp_u=length(CP{ip}(:,1,1));
        cp_v=length(CP{ip}(1,:,1));
    
        % read solution file
        fname=strcat(loc,'p',sprintf('%02d',ip),'_kin/t',sprintf('%06d',n),'.h5');
        [CP_pot,CP_wrec,CP_ii,CP_ia]=getdatah5(fname,cp_u,cp_v);

        % save variable for action potential/ionic current plot
        if (ip==2)
            vv(n+1) = get_point_eval(p(ip),0.5,U{ip},q(ip),0.5,V{ip},CP{ip},CP_pot);
            iion(n+1) = get_point_eval(p(ip),0.5,U{ip},q(ip),0.5,V{ip},CP{ip},CP_ii);

        end        
       
        % evaluate potential on a subgrid to get a smooth scalar field
        switch var
            case 1
                C = create_cont(p(ip),q(ip),U{ip},V{ip},CP{ip},CP_pot,nn);                
                titi='action potential v [mV]';                
            case 2
                C = create_cont(p(ip),q(ip),U{ip},V{ip},CP{ip},CP_wrec(:,:,1),nn);
                titi='1st gating variable';
            case 3
                C = create_cont(p(ip),q(ip),U{ip},V{ip},CP{ip},CP_ii,nn);
                titi='ionic current [mA/cm$^2$]';
            case 4
                C = create_cont(p(ip),q(ip),U{ip},V{ip},CP{ip},CP_ia,nn);
                titi='applied current [mA/cm$^2$]';
        end

        nexttile(1)
        surf(X(:,:,ip),Y(:,:,ip),Z(:,:,ip),C,'EdgeColor','none')
        hold on

    end

    PP = get_point_surf(p(2),0,0.5,U{2},q(2),0,0.5,V{2},CP{2});

    delete(findall(gcf,'type','annotation'))

    % ventricle plot
    colorbar
    clim([-80,20])    
    axis equal
    axis([-4,4,-4,4,-7,3])
    scatter3(PP(1),PP(2),PP(3),20,'MarkerEdgeColor','none','MarkerFaceColor','r')
    xlabel('x [cm]','FontSize',14,'Interpreter','latex');
    ylabel('y [cm]','FontSize',14,'Interpreter','latex');
    zlabel('z [cm]','FontSize',14,'Interpreter','latex');
    tit=strcat(titi,', t =',sprintf('%4.2f',n*dtfig),' [ms]'); 
    title(tit,'FontSize',14,'Interpreter','latex')
    camlight
    lightangle(100,30)
    view(150,20)
    hold off  

    % action potential plot
    nexttile(2)
    yyaxis left
    plot(linspace(0,t(n+1),length(vv)),vv,'-')
    hold on
    axis([0,550,-90,30])
    ylabel('v$_P$ [mv]','FontSize',14,'Interpreter','latex');
    xlabel('t [ms]','FontSize',14,'Interpreter','latex');
    yyaxis right
    plot(linspace(0,t(n+1),length(vv)),iion,'-')   
    ylim([-1,8].*10^-5)
    ylabel('(i$_{ion})_P$ [mA/cm$^2$]','FontSize',14,'Interpreter','latex');
    hold off  

    annotation(gcf,'textbox',[0.339,0.548,0.028,0.060],...
    'Color','r','String',{'P'},'LineStyle','none','Interpreter','latex',...
    'FontSize',16,'FitBoxToText','off');

    drawnow    
    disp(['step ' num2str(n)])
    

    if (wrf)
        picname=strcat('./frames/frame_',sprintf('%05d',n));
        print('-dpng','-r300',picname);                        
    end    
    
end

