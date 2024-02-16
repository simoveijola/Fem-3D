%% Best Placement of Holes in Cylinder to Achieve Target Surface Temperature - A Parametric Study
% Copyright (c) 2015, MathWorks, Inc.
%
% This examples conducts a parametric study in which heat conduction
% simulation is performed over a set of similar geometries to determine
% which geometry "best" meets an average nodal target temperature on a
% specified output face. The geometry is a cylinder like structure as shown
% below and has a ring of holes running longitudinally through the
% structure. The problem has the following characteristics:
%
% 
% <<annotated_geometry.png>>
% 
% *Boundary conditions*
% 
% * The input heat source is applied on the faces of the holes
% * The longitudinal surface (called 'output face' later in text) and the
% surface on the center protrusion have convective boundary conditions. The
% objective of the example is to achieve a target average nodal temperature on
% output face in the "best" manner possible. 
% * All other faces not indicated above are insulated and thus have zero
% Neumann boundary conditions.
% 
% *Geometry*
% 
% * Each geometry has a unique pair of (#holes, radius of ring of holes).
% All other geometry parameters are held constant.
% * Although it is possible to exploit symmetry in some of the geometries
% in order to reduce the problem to 2 geometry dimensions, it was not done
% so in this example.
% 
% *Results*
% 
% * Results are collected from all simulations and the best geometry in
% terms of lowest max-min temperature spread on the output face and the
% best geometry for lowest operating cost (input flux) are identified.
% * The implementation uses 'parfeval' from the Parallel Computing Toolbox
% to speed up the parametric study.
% 

function efficient_heating

%% Import candidate geometries and create models
% The STL files are read - each file corresponds to a different
% parameter-pair: (#holes, radius of ring of holes).
fileList = ls('cyl_*.STL');
fileList = mat2cell(fileList,ones(size(fileList,1),1));
%%
% The PDE is a scalar, laplace equation
N = 1;
%%
% A table will be created to organize and report results from all the runs.
Model = cellfun(@(~) createpde(N),fileList,'UniformOutput',false);
% Table columns corresponding to #holes and radius (note: this is the radius of
% the ring and *not* the radius of the hole) are extracted.
paramList = cellfun(@(fileName) regexpi(fileName,'cyl_(.*)_(.*).STL','tokens'),fileList, 'UniformOutput',false);
NumHoles = cellfun(@(entry) str2double(entry{1}(1)),paramList,'UniformOutput',false);
HolesRadius = cellfun(@(entry) str2double(entry{1}(2)),paramList,'UniformOutput',false);
%%
% Create table
T = [table(Model), cell2table(NumHoles,'RowNames',fileList), cell2table(HolesRadius)];
% Create cell array of geometries. Each geometry gets its own model
for k = 1:size(T,1)
    importGeometry(T.Model{k},T.Properties.RowNames{k});
    % The relation of faces to holes is known; report errors for unexpected relation 	
    if T.Model{k}.Geometry.NumFaces ~= (3 + T.NumHoles(k) + 2)
        error('unexpected number of faces');
    end
end
%%
% Sorted table, first by #holes and then by radius
T = sortrows(T,{'NumHoles','HolesRadius'},{'ascend','ascend'});
%%
% Plot two extreme geometries to show the range of geometry variations
figure
pdegplot(T.Model{1},'FaceLabels','on');
title(T.Properties.RowNames{1});
view(0,90);
figure
pdegplot(T.Model{end},'FaceLabels','on');
title(T.Properties.RowNames{end});
view(0,90);

%% Input and Output faces for all geometries
%%
% Any face in |(inputFacesBegin:(inputFacesBegin +numHoles))| is an input
% heat source face
inputFacesBegin = 4;
%%
% Output face on which the average nodal temperature is measured. Note that
% the output is *not* a true average temperature and is only an average of
% nodes on the output face and is therefore more mesh-dependent than we
% would like and is an approximation. This approximation is used as a proxy
% for the true output to keep the example relatively simple.
outputFace = 1;

%% Solution and desired output setup
%%
% Table column for operating cost (total flux going into solid via the
% input heat source faces); it is desirable to minimize this
OperatingCost = zeros(size(T,1),1);
%%
% |Constant| and |MaxMinSpread| columns correspond to solutions for the
% corresponding constant and variables contributions of the affine boundary
% conditions. Solving for the solution this way i.e. two times per
% geometry, will let us scale the variable part to match the target average
% temperature.

MaxMinSpread = zeros(size(T,1),1);
%%
% Table column for scale factor for |MaxMinSpread|
ScaleForTargetMaxMinSpread = zeros(size(T,1),1);
%%
% Add aforementioned additional columns to table
T = [T table(MaxMinSpread,ScaleForTargetMaxMinSpread,OperatingCost)];

%% Input (non-geometry) setup
%%
% Ambient temperature
ambientTemp = 6;
%%
% Target average nodal temperature on output face
targetMaxMinSpread = 3;
%%
% PDE coefficients for laplace equation (heat conduction)
c = 1e-1;
a = 0;
f = 0;
%%
% Boundary conditions related input
convectiveHeatTransferCoeff = 0.3;

%% Function for applying boundary conditions, meshing, and solving per geometry
    function [maxMinSpread,resultVariableBC] = solveGeometry(model,numHoles)
        % generate mesh with 'hmax' = 1/4th of hole radius
        model.generateMesh('hmax',0.25/4);
        % extract nodes on output face
        [~,e,~] = meshToPet(model.Mesh);
        outputFaceNodes = e.getNodes(outputFace);
        % variable component of boundary conditions
        % generalized Neumann BC on output face and also face on center
        % protrusion
        model.applyBoundaryCondition('Face',[outputFace,model.Geometry.NumFaces],...
            'q',convectiveHeatTransferCoeff);
        % apply unit flux on input heat source faces
        model.applyBoundaryCondition('Face',(inputFacesBegin:(inputFacesBegin + numHoles)),'g',1);
        resultVariableBC = assempde(model,c,a,f);
        temp = resultVariableBC(outputFaceNodes);
        maxMinSpread = max(temp) - min(temp);
    end

%% Solve for all geometries
%%
% Use parfeval to perform asynchronous computation and speed up overall
% simulation time
pool = gcp();
%%
% Futures are created for the geometries
for idx = 1:size(T,1)
    % coefficients of PDE
    F(idx) = parfeval(pool,@solveGeometry,2,T.Model{idx},T.NumHoles(idx));
end
%%
% Populate table with results of computations that are performed
% asynchronously
for idx = 1:size(T,1)
    % get result for next geometry that was solved
    [completedIdx,maxMinSpread] = fetchNext(F);
    % compute average nodal temperatures for the variable and constant
    % components of the output. As was been noted previously, this is just
    % a proxy for the true average temperature on the output face.
    T.MaxMinSpread(completedIdx) = maxMinSpread;
    fprintf('%d of %d models simulated\n',completedIdx,size(T,1));
    % compute scale factor for variable output.
    T.ScaleForTargetMaxMinSpread(completedIdx) = targetMaxMinSpread./T.MaxMinSpread(completedIdx);
end

%% Report and visualize results
%%
T.OperatingCost = T.ScaleForTargetMaxMinSpread.*T.NumHoles;
%%
% Top-5 operating cost sorted from smallest to largest
TOpCost = sortrows(T,'OperatingCost');
TOpCost(1:5,:)
[~,resultVariableBC] = solveGeometry(TOpCost.Model{1},TOpCost.NumHoles(1));
u_Optimal = resultVariableBC*TOpCost.ScaleForTargetMaxMinSpread(1) + ambientTemp;
figure
pdeplot3D(TOpCost.Model{1},'colormapdata',u_Optimal);
view(45,90);
snapnow
view(114,51);
snapnow
end