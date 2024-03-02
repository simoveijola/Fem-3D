%% quad on unit tetrahedron
w = 1/24*ones(1,4);
x1=[0.1381966011250105,0.1381966011250105,0.1381966011250105]; x2=[0.5854101966249685,0.1381966011250105,0.1381966011250105];
x3=[0.1381966011250105,0.5854101966249685,0.1381966011250105]; x4=[0.1381966011250105,0.1381966011250105,0.5854101966249685];
ti = [x1;x2;x3;x4]';

%% affine mapping from unit tetrahedron to arb tetrahedron: see below


%% contributions to Ahat on Ti
% For each tetrahedron, we have four basis functions with support 
% on this tetrahedron. These linear basis functions correspond to the nodes
% of T_i, and the indeces are given by the connectivity matrix of the 3D
% discretization. 

%% Volume of T_i: see below


%% Reference basis function on unit tetrahedron:
% we assume that the origin is mapped to the first node in Ti
phi(:,1) = 1 - ti(1,:) - ti(2,:) - ti(3,:);
phi(:,2) = ti(1,:);
phi(:,3) = ti(2,:);
phi(:,4) = ti(3,:);
dphi = [-1 -1 -1 ;
        1 0 0 ;
        0 1 0 ;
        0 0 1];
% row i is gradient of phi(:,1)

%% funcs
function a = AT(nodes)
    n1 = nodes(:,1);
    M = [nodes(:,2)-n1, nodes(:,3)-n1, nodes(:,4)-n1];
    b = n1*ones(nodes);
    a = M*nodes + b;
end

function V = volume(nodes)
    n4 = nodes(:,4);
    M = [nodes(:,1)-n4, nodes(:,2)-n4, nodes(:,3)-n4];
    V = 1/6*abs(det(M));
end
