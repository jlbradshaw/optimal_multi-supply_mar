% Copyright Jonathan L. Bradshaw 2018

% function pipelines = pipeline_routine(wrf_array,sbasin_array)
function [pipeLength_ft, pipeFeasible_index, pump_reqs, prv_reqs] = pipeline_routine(elev_profile_m, plantCapacity_afy, po)

% Purpose: Determine which pipe size options are feasible and compute the
% pumping requirements.
% Input: elevation profile of the pipeline in m, plant capacity in AFY,
% pipe options (po)
% Output: pipeline length in feet, diameter and Hazen William coefficients of feasible pipes,
% required pumps for each feasible pipe diameter, required pressure reducing values (PRVs)

% Maximum and minimum pressure head allowed in the pipeline
psiMax_ft = 461.33; % 200 pounds per square inch of pressure as head in feet
psiMin_ft = 461.33/20.; % 10 pounds per square inch of pressure as head in feet

% Min and max permitted peak pipe velocity
velMin_fps = 1; % minimum permitted peak velocity in the pipe in feet per second
velMax_fps = 7; % maximum permitted peak velocity in the pipe in feet per second

d = po.diameters; % assume diameters is in units of FEET
% For each pipe diameter, compute the min and max permitted flowrate
area = pi()*(d/2).^2; 
flowMin_afy = units_check(velMin_fps*area,'cfs2afy');
flowMax_afy = units_check(velMax_fps*area,'cfs2afy');

% Determine which pipe diameters and roughnesses are within the feasible
% flow range based on the peak flow rate (i.e., plant Capacity)
pipeFeasible_index = (plantCapacity_afy >= flowMin_afy) & (plantCapacity_afy <= flowMax_afy);

diam_feasible_ft = po.diameters(pipeFeasible_index);
C_feasible = po.Cs(pipeFeasible_index);

% Compute a vector of the head loss per unit pipe length due to friction
% for the plant capacity and feasible diameters
phi = 4.73; % assumes US units
hfs_perL_ft = 1.1*(phi.*units_check(plantCapacity_afy,'afy2cfs')^1.85)./((C_feasible.^1.85).*(diam_feasible_ft.^4.87)); % compute the head loss due to fricition per unit length (here, feet) of pipeline (i.e., L = 1). Each row is a different pipe type (i.e., C and d), and each column is a different flow rate. The coefficient of 1.1 is to account for minor losses, which is assumed to be 10% of major frictional losses.

% Create a cell array that contains the pumping requirements for each of the
% feasible pipe options
pump_reqs = cell(length(diam_feasible_ft),1);
prv_reqs = pump_reqs;
for i = 1:length(diam_feasible_ft)
    [pipeLength_ft, pump_reqs{i}, prv_reqs{i}] = pumpdesign(elev_profile_m, hfs_perL_ft(i), psiMax_ft, psiMin_ft);
end