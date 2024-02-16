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
% <<publish_pic.png>>
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

function [TOpCost,TMaxMinSpread] = heating_simplified

%% Import candidate geometries and create models
% The STL files are read - each file corresponds to a different
% parameter-pair: (#holes, radius of ring of holes).

fileList = ls('cyl_*.STL');
fileList = mat2cell(fileList,ones(size(fileList,1),1));
fileList = fileList(1:2);
%%
% The PDE is a scalar, laplace equation
N = 1;

%%
% A table will be created to organize and report results from all the runs.
% Columns corresponding to #holes and radius (note: this is the radius of
% the ring and *not* the radius of the hole) are extracted.
paramList = cellfun(@(fileName) regexpi(fileName,'cyl_(.*)_(.*).STL','tokens'),fileList, 'UniformOutput',false);
NumHoles = cellfun(@(entry) str2double(entry{1}(1)),paramList,'UniformOutput',false);
HolesRadius = cellfun(@(entry) str2double(entry{1}(2)),paramList,'UniformOutput',false);
%%
% Create table sorted by #holes and radius and in that order
T = [cell2table(NumHoles,'RowNames',fileList), cell2table(HolesRadius)];
T = sortrows(T,{'NumHoles','HolesRadius'},{'ascend','ascend'});
%%
% Column corresponding to the indices of the files
T.GeometryIndex = (1:length(fileList))';
%%
% Give |fileList| identical ordering to that of the table
fileList = fileList(T.GeometryIndex);
%%
% Create cell array of geometries. Each geometry gets its own model
modelList = cellfun(@(~) createpde(N),fileList,'UniformOutput',false);
for k = 1:length(modelList)
    importGeometry(modelList{k},fileList{k});
    % The relation of faces to holes is known; report errors for unexpected relation 	
    if modelList{k}.Geometry.NumFaces ~= (3 + T.NumHoles(k) + 2)
        error('unexpected number of faces');
    end
end
%%
% Plot two extreme geometries to show the range of geometry variations
figure
pdegplot(modelList{1},'FaceLabels','on');
title(fileList{1});
view(0,90);
figure
pdegplot(modelList{end},'FaceLabels','on');
title(fileList{end});
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
%
%%
% Table columns for capturing max, min temperatures and spread on output
% face; it is desirable to have a low spread
MinTemp = zeros(length(modelList),1);
MaxTemp = zeros(length(modelList),1);
MaxMinSpread = zeros(length(modelList),1);
%%
% Table column for operating cost (total flux going into solid via the
% input heat source faces); it is desirable to minimize this
OperatingCost = zeros(length(modelList),1);
%%
% |Constant| and |Variable| columns correspond to solutions for the
% corresponding constant and variables contributions of the affine boundary
% conditions. Solving for the solution this way i.e. two times per
% geometry, will let us scale the variable part to match the target average
% temperature.
Constant = zeros(length(modelList),1);
Variable = zeros(length(modelList),1);
%%
% Table column for scale factor for |Variable|
InputForTargetTemp = zeros(length(modelList),1);
%%
% Add aforementioned additional columns to table
T = [T table(Constant,Variable,InputForTargetTemp,MinTemp,MaxTemp,MaxMinSpread,OperatingCost)];
%%
% |u| will hold all the solutions
u = cell(length(modelList),1);

%% Input (non-geometry) setup
%%
% Ambient temperature
ambientTemp = 6;
%%
% Target average nodal temperature on output face
targetTemp = 15;
%%
% PDE coefficients for laplace equation (heat conduction)
c = 1e-1;
a = 0;
f = 0;
%%
% Boundary conditions related input
convectiveHeatTransferCoeff = 0.3;

%% Function for applying boundary conditions, meshing, and solving per geometry
    function [resultVariableBC,resultConstantBC,outputFaceNodes] = solveGeometry(model,numHoles)
        % generate mesh with 'hmax' 1/4th of hole radius
        model.generateMesh('hmax',0.25/4);
        % extract nodes on output face
        [~,e,~] = meshToPet(model.Mesh);
        outputFaceNodes = e.getNodes(outputFace);
        % variable component of boundary conditions
        % generalized Neumann BC on output face and also face on center
        % protrusion
        model.applyBoundaryCondition('Face',[outputFace,modelList{idx}.Geometry.NumFaces],...
            'q',convectiveHeatTransferCoeff);
        % apply unit flux on input heat source faces
        model.applyBoundaryCondition('Face',(inputFacesBegin:(inputFacesBegin + numHoles)),'g',1);
        % *** resultVariableBC = assempde(model,c,a,f);        
        resultVariableBC = numHoles;
        % constant component of boundary conditions is calculated by
        % subtracting the variable component calculated above from the
        % output for the full set of boundary conditions
        model.applyBoundaryCondition('Face',[outputFace,model.Geometry.NumFaces],...
            'q',convectiveHeatTransferCoeff,'g',convectiveHeatTransferCoeff*ambientTemp);
        % *** resultConstantBC = assempde(model,c,a,f) - resultVariableBC;  
        resultConstantBC = 2*numHoles;
                
    end
%% Solve for all geometries
%%
% Use parfeval to perform asynchronous computation and speed up overall
% simulation time
pool = gcp();
%%
% Futures are created for each geometry
for idx = 1:length(modelList)
    % coefficients of PDE
    F(idx) = parfeval(pool,@solveGeometry,3,modelList{idx},T.NumHoles(idx));
end
%%
% Populate table with results of computations that are performed
% asynchronously
for idx = 1:length(modelList)
    % get result for next geometry that was solved
    [completedIdx,resultVariableBC,resultConstantBC,outputFaceNodes] = fetchNext(F);
    % compute average nodal temperatures for the variable and constant
    % components of the output. As was been noted previously, this is just
    % a proxy for the true average temperature on the output face.
%{    
    % *** 
    T.Variable(completedIdx) = mean(resultVariableBC(outputFaceNodes));
    T.Constant(completedIdx) = mean(resultConstantBC(outputFaceNodes));
    fprintf('%d of %d models simulated\n',completedIdx,length(modelList));
    % compute scale factor for variable output.
    T.InputForTargetTemp(completedIdx) = (targetTemp-T.Constant(completedIdx))./T.Variable(completedIdx);
    % compute actual result
    u{completedIdx} = resultVariableBC*T.InputForTargetTemp(completedIdx) + resultConstantBC;
    % compute min max temperatures on output face
    T.MinTemp(completedIdx) = min(u{completedIdx}(outputFaceNodes));
    T.MaxTemp(completedIdx) = max(u{completedIdx}(outputFaceNodes));
    %}
end

%% Report and visualize results
%%
% Calculate max-min spread and operating cost
T.MaxMinSpread = T.MaxTemp - T.MinTemp;
T.OperatingCost = T.InputForTargetTemp.*T.NumHoles;
%%
% Top-5 operating cost sorted from smallest to largest
TOpCost = sortrows(T,'OperatingCost');
TOpCost(1:2,:)
%%
% The operating cost numbers are all close to each other so it isn't
% effective as a criterion.
%%
% Top-5 max-min spread sorted from smallest to largest
TMaxMinSpread = sortrows(T,'MaxMinSpread');
TMaxMinSpread(1:2,:)
%%
% So, numerous small radii holes instead of fewer large radii holes and the
% holes that are closest to the center affect the temperature in the
% direction we want. The solution with the smallest max-min spread on the
% output face and one that matches the target average temperature is
% plotted below.
figure
%*** pdeplot3D(modelList{TMaxMinSpread.GeometryIndex(1)},'colormapdata',u{TMaxMinSpread.GeometryIndex(1)});
view(45,90);
end