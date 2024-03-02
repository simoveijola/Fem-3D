%% how to import geometry and initialize mesh:
model = createpde;
importGeometry(model, "3D_model_of_a_Cube.stl");
generateMesh(model, "GeometricOrder","linear")

p = model.Mesh.Nodes; t = model.Mesh.Elements;
% from the mesh we get the nodes and the connectivity matrix.
pdeplot3D(model,'NodeLabels','off', 'FaceAlpha',0.3);

%% how to plot:
pdeplot3D(model.Mesh, ColorMapData=ones(size(p,2),1))

%% cube:
scatter3(p(1,:),p(2,:),p(3,:))
idof = min(p,[],1)==0 | max(p,[],1)==30;
sum(idof)
extnodes = p(:,idof);
scatter3(extnodes(1,:),extnodes(2,:),extnodes(3,:))

%%
outer_nodes = solve_outer_nodes(t);
%%

scatter3(p(1,outer_nodes),p(2,outer_nodes),p(3,outer_nodes))
%%
f = @(x) (1);
idof = setdiff(t, outer_nodes);
w = @(x) 15^(-6)*x(1)*(30-x(1))*x(2)*(30-x(2))*x(3)*(30-x(3)) - 20;
dw{1} = @(x) 15^(-6)*(30-2*x(1))*x(2)*(30-x(2))*x(3)*(30-x(3));
dw{2} = @(x) 15^(-6)*x(1)*(30-x(1))*(30-2*x(2))*x(3)*(30-x(3));
dw{3} = @(x) 15^(-6)*x(1)*(30-x(1))*x(2)*(30-x(2))*(30-2*x(3));
uh = threedsolver(p,t,f,[1,0,0;0,1,0;0,0,1], idof,w,dw);

plot3D_solution(p,t,uh)