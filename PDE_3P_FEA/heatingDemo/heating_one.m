%% Step 1. Import Geometry
% Copyright (c) 2015, MathWorks, Inc.
myPde = createpde(1); 
importGeometry(myPde,'cyl_7_1.0.STL');
figure
pdegplot(myPde,'FaceLabels','on');

%% Step 2. Specify PDE coefficients
% PDE coefficients for laplace equation (heat conduction)
c = 1e-1;
a = 0;
f = 0;

%% Step 3. Specify BC
% Ambient temperature
ambientTemp = 6;
% Convective heat transfer coefficient
hc = 0.3;
% Faces with convective boundary conditions
sideFace = 1;
protrusionFace = myPde.Geometry.NumFaces;
% Apply convective boundary condition
myPde.applyBoundaryCondition('Face',[sideFace,protrusionFace],...
    'q',hc,'g',ambientTemp);
% Faces with heat source
inputFaces = (4:11);
% Apply unit flux on input heat source faces
myPde.applyBoundaryCondition('Face',inputFaces,'g',1);

%% Step 4. Mesh
% Generate mesh with 'hmax' = 1/4th of hole radius
myPde.generateMesh('hmax',0.25/2);

%% Solve
result = assempde(myPde,c,a,f);

%% Visualize
figure
pdeplot3D(myPde,'colormapdata',result);