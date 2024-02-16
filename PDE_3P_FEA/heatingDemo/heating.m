%% Best Placement of Holes in Cylinder to Achieve Target Average Temperature - A Parametric Study
% Copyright (c) 2015, MathWorks, Inc.
%
% This examples conducts a parametric study in which heat conduction
% simulation is performed over a set of similar geometries to determine
% which geometry "best" meets an average temperature on an specified output
% area. The geometry is a cylinder like structure as shown below and has a
% ring of holes running longitudinally through the structure. The problem
% has the following characteristics:
%
% 
% <<publish_pic.png>>
% 
% *Boundary conditions*
% 
% * The input heat source is applied on the faces of the holes.
% * The longitudinal surface and the surface on the center protrusion have
% convective boundary conditions. The output surface is the rectangular
% subsurface of the longitudinal surface. The objective of the example is
% to achieve a target average temperature on output surface in the "best"
% manner possible as described later.
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
% terms of lowest max-min temperature spread on the longitudinal face and
% the best geometry for lowest operating cost (input flux) are identified.
% * The implementation uses 'parfeval' from the Parallel Computing Toolbox
% to speed up the parametric study.
% 
function heating

%% Import candidate geometries and create models
%%
% The STL files are read - each file corresponds to a different
% parameter-pair: (#holes, radius of ring of holes).
fileList = ls('cyl_*.STL');
fileList = mat2cell(fileList,ones(size(fileList,1),1));
%%
% The PDE is a scalar, laplace equation
N = 1;
%%
% A table will be created to organize data and results from all runs
%%
% Table column corresponding to the PDE models
Model = cellfun(@(~) createpde(N),fileList,'UniformOutput',false);
%%
% Table columns corresponding to #holes and radius (note: this is the
% radius of the ring and *not* the radius of the hole) are extracted.
paramList = cellfun(@(fileName) regexpi(fileName,'cyl_(.*)_(.*).STL','tokens'),fileList, 'UniformOutput',false);
NumHoles = cellfun(@(entry) str2double(entry{1}(1)),paramList,'UniformOutput',false);
HolesRadius = cellfun(@(entry) str2double(entry{1}(2)),paramList,'UniformOutput',false);
%%
% Create table
T = [table(Model), cell2table(NumHoles,'RowNames',fileList), cell2table(HolesRadius)];
%%
% Import geometries into the PDE models
for k = 1:size(T,1)
    importGeometry(T.Model{k},T.Properties.RowNames{k});
    % The relation of faces to holes is known; report errors for unexpected
    % relation
    if T.Model{k}.Geometry.NumFaces ~= (3 + T.NumHoles(k) + 2)
        error('unexpected number of faces');
    end
end
%%
% Sort table, first by #holes and then by radius
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

%% Input setup
%%
% Ambient temperature
ambientTemp = 6;
%%
% Target average nodal temperature on output surface
targetTemp = 15;
%%
% PDE coefficients for laplace equation (heat conduction)
c = 1e-1;
a = 0;
f = 0;
%%
% Any face in |(inputFacesBegin:(inputFacesBegin +numHoles))| is an input
% heat source face
inputFacesBegin = 4;

%% Output setup
%%
% Table column for capturing max-min temperature on output
% area; it is desirable to have a low spread
MaxMinSpread = zeros(size(T,1),1);
%%
% Table column for operating cost (total flux going into solid via the
% input heat source faces); it is desirable to minimize this
OperatingCost = zeros(size(T,1),1);
%%
% |AvgTempVariable| column corresponds to the variable contribution towards
% the average temperature solution on the output surface. The constant
% contribution is simply |ambientTemp|.
AvgTempVariable = zeros(size(T,1),1);
%%
% Table column for scale factor for |AvgTempVariable| to help match |targetTemp|
InputForTargetTemp = zeros(size(T,1),1);
%%
% Add these columns to table
T = [T table(AvgTempVariable,InputForTargetTemp,MaxMinSpread,OperatingCost)];
%%
% Face on which the max-min temperature spread is measured. 
MaxMinSpreadFace = 1;
%%
% Output surface in XZ plane: |-offsetX:offsetX, offsetY, minZ:maxZ| where average
% temperature is calculated
offsetY = -1.875;
offsetX = sqrt(2^2-offsetY^2);
minZ = 0;
maxZ = 1;
%%
% Convective heat transfer coefficient
hc = 0.3;

%% Function for applying boundary conditions, meshing, and solving per geometry
    function [avgTemp,maxMinSpread,resultVariableBC] = solveGeometry(model,numHoles)
        % generate mesh with 'hmax' = 1/4th of hole radius
        model.generateMesh('hmax',0.25/4);
        % extract nodes on MaxMinSpreadFace
        [p,e,t] = meshToPet(model.Mesh);
        maxMinSpreadFaceNodes = e.getNodes(MaxMinSpreadFace);
        % variable component of boundary conditions is generalized Neumann BC
        % on maxMinSpreadFaceNodes and also face on center protrusion
        model.applyBoundaryCondition('Face',[MaxMinSpreadFace,model.Geometry.NumFaces],...
            'q',hc);
        % apply unit flux on input heat source faces
        model.applyBoundaryCondition('Face',(inputFacesBegin:(inputFacesBegin + numHoles)),'g',1);
        % solve to get result for variable BC
        resultVariableBC = assempde(model,c,a,f);
        % calculate max-min spread
        t1 = resultVariableBC(maxMinSpreadFaceNodes);
        maxMinSpread = max(t1) - min(t1);
        % calculate average temp. on rectangular output surface area
        myInterpolant = pdeInterpolant(p,t,resultVariableBC);
        function res = intFun(x,z)
            % Y offset is pushed a bit inwards to avoid missing data along
            % Y axis
            res = evaluate(myInterpolant,x,(offsetY+0.01)*ones(size(x)),z);
            res = reshape(res,size(x));
            % NaNs dues to XZ plane overshoot are set to zero
            res(find(isnan(res)))=0;
            return
        end
        area = 2*offsetX*(maxZ-minZ);
        avgTemp = integral2(@intFun,-offsetX,+offsetX,minZ,maxZ)/area;
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
    F(idx) = parfeval(pool,@solveGeometry,3,T.Model{idx},T.NumHoles(idx));
end
%%
% Populate table with results of computations that are performed
% asynchronously
for idx = 1:size(T,1)
    % get result for next geometry that was solved
    [completedIdx,avgTemp,maxMinSpread] = fetchNext(F);
    % Set the average temperature contribution of the variable part of the
    % BC
    T.AvgTempVariable(completedIdx) = avgTemp;
    fprintf('%d of %d models simulated\n',completedIdx,size(T,1));
    % compute scale factor for the contribution of the variable part. As
    % mentioned earlier ambientTemp corresponds to contribution of the constant part.
    T.InputForTargetTemp(completedIdx) = (targetTemp-ambientTemp)./T.AvgTempVariable(completedIdx);
    T.MaxMinSpread(completedIdx) = maxMinSpread;
end

%% Report and visualize results
%%
% Calculate operating cost
T.OperatingCost = T.InputForTargetTemp.*T.NumHoles;
%%
% Top-5 operating cost sorted from smallest to largest
TOpCost = sortrows(T,'OperatingCost');
TOpCost(1:5,:)
%%
% Plot result for geometry with lowest operating cost
[~,~,resultVariableBC] = solveGeometry(TOpCost.Model{1},TOpCost.NumHoles(1));
u_Optimal = resultVariableBC*TOpCost.InputForTargetTemp(1) + ambientTemp;
figure
pdeplot3D(TOpCost.Model{1},'colormapdata',u_Optimal);
view(45,90);
snapnow
view(114,51);
snapnow
%%
% Top-5 max-min spread sorted from smallest to largest
TMaxMinSpread = sortrows(T,'MaxMinSpread');
TMaxMinSpread(1:5,:)
%%
% Plot result for geometry with lowest max-min spread
[~,~,resultVariableBC] = solveGeometry(TMaxMinSpread.Model{1},TMaxMinSpread.NumHoles(1));
u_Optimal = resultVariableBC*TMaxMinSpread.InputForTargetTemp(1) + ambientTemp;
figure
pdeplot3D(TMaxMinSpread.Model{1},'colormapdata',u_Optimal);
view(45,90);
snapnow
view(114,51);
snapnow

%% Takeaways
%%
% 
% * Programmatically solve PDEs and perform parametric studies. Useful if
% there is little fundamental variation between different design points.
% * Perform custom post-processing with useful reports
% * Use Parallel Computing Toolbox to accelerate simulations
% 

end