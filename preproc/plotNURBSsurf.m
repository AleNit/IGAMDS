function plotNURBSsurf(p,q,U,V,CP,el,conp)
% plots the surface, elements and control points
% J. Kiendl

[X,Y,Z] = create_surf(p,q,U,V,CP);

% geometry
surf(X,Y,Z,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
hold on

if (el) % element edges
    create_el_edges(p,q,U,V,CP)
end

if (conp)   % control points and polygon
    create_conpolygon(CP)
end

% camlight left; lighting phong;
axis equal;
xlabel('x [cm]','FontSize',14,'Interpreter','latex');
ylabel('y [cm]','FontSize',14,'Interpreter','latex');
zlabel('z [cm]','FontSize',14,'Interpreter','latex');

end
