% Copyright Jonathan L. Bradshaw 2018

function [A, b, Aeq, beq, lb, ub, rw_indices, sw_indices, infiltration_indices, storage_indices, x0] = getConstraints_papertwo_twoWRFseries(rwCapacity_mgd_point, number_timesteps, wasteFractionMFRO, ET_mgd, wwAvailable_mgd, swAvailable_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, storageCapacity_mg, infilCapacity_mgd)

numWRFs = 2;

% if we are not including ET in constraints, set ET to 0
if ~is_ETinConstraints
    ET_mgd = ET_mgd*0;
end


%% Indices
number_processes = 0;                                                        % total number of processes represented

% Flows
number_flows = 2;
flow_cell = cell(1,number_flows);           % Create a 1-D array to store the flow indices

% loop through number of flows and create indices to assign to the
% different flows.
for i = 1:number_flows
    flow_cell{i} = (1:1:number_timesteps)+(i-1)*number_timesteps + number_processes;
end
rw_indices1 = flow_cell{1};
rw_indices2 = flow_cell{2};
rw_indices = cell(1,2);
rw_indices{1} = rw_indices1;
rw_indices{2} = rw_indices2;

temp = number_flows*number_timesteps+number_processes+1;
number_dummy_var = 0;

%SW usage
sw_indices=(temp:1:temp+number_timesteps-1);

%Infiltration rate
infiltration_indices=(sw_indices(end)+1:1:sw_indices(end)+number_timesteps);

% Storage
storage_indices = (infiltration_indices(end)+1:1:infiltration_indices(end)+number_timesteps);

number_decisionVars = number_processes + number_timesteps*number_flows + number_dummy_var + length(sw_indices) + length(infiltration_indices) + length(storage_indices); % the number of decision variables


%% Building the matrix for inequalities

A_timestepTemplate = sparse(zeros(number_timesteps, number_decisionVars));                    %Ttemplate for contraint matrices.

ind_rw1_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,rw_indices1);                    %Creates linear index for q0
ind_rw2_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,rw_indices2);                    %Creates linear index for q0
ind_sw_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,sw_indices);
ind_infil_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,infiltration_indices);
ind_storage_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,storage_indices);


%% Constraints
Aeq1 = A_timestepTemplate;

% Special case for the first time step: assume no initial storage
for i = 2:length(storage_indices) % for loop starts at i = 2 to avoid reference to storage(t=0), which doesn't exist
    Aeq1(i, storage_indices(i-1)) = 1;
end

Aeq1(ind_sw_A) = 1;
Aeq1(ind_rw1_A)= 1;
Aeq1(ind_rw2_A)= 1;

Aeq1(ind_infil_A) = -1;
Aeq1(ind_storage_A) = -1;

beq1 = ET_mgd;

%% Iterate through RW capacities to complete constraints and run optimization process


%% Bounds

lb_rw1_mgd = zeros(1,number_timesteps);
lb_rw2_mgd = zeros(1,number_timesteps);

lb_sw_mgd = zeros(1,length(sw_indices));
lb_infil_mgd = zeros(1,length(infiltration_indices));
lb_storage_mg = zeros(1,length(storage_indices));

ub_storage_mg = storageCapacity_mg*ones(1,number_timesteps);
ub_infil_mgd = infilCapacity_mgd*ones(1,number_timesteps);
ub_swWIntake_mgd = min(swAvailable_mgd',intakeLimits_mgd);

A = [];
b = [];

Aeq = sparse(Aeq1);
beq = sparse(beq1);

ub_rw1_mgd = min(wwAvailable_mgd{1}, rwCapacity_mgd_point(1));
ub_rw2_mgd = min(wwAvailable_mgd{2}, rwCapacity_mgd_point(2));

if is_uncap_fixed % if the unused capacity is supposed to be fixed (e.g., by month)
    
    % replace the upper bound of the rw input with the alternative unused capacity
    ub_rw1_mgd = min(uncap_daily_mgd, ub_rw1_mgd);
    ub_rw2_mgd = min(uncap_daily_mgd, ub_rw2_mgd);
    if ~is_minRWCapacity_fixed
        disp('Need to complete for multiple WRFs')        
    end
    
    
end

% Compile the lower and upper bounds
lb = sparse([lb_rw1_mgd, lb_rw2_mgd, lb_sw_mgd, lb_infil_mgd, lb_storage_mg]);
ub = sparse([ub_rw1_mgd, ub_rw2_mgd, ub_swWIntake_mgd, ub_infil_mgd, ub_storage_mg]);

%% Compute an initial guess for interior point method
swInput0 = zeros(1,number_timesteps);
storage0 = zeros(1,number_timesteps);
rwProduction1_0 = lb_rw1_mgd;
rwProduction2_0 = lb_rw2_mgd;
rwProduction0_cell = {rwProduction1_0, rwProduction2_0};
infiltration0 = zeros(1,number_timesteps);

initial_storage = 0;

for i = 1:number_timesteps
    if i == 1
        temp_storage = initial_storage;
    else
        temp_storage = storage0(i-1);
    end
    
    swInput0(i) = min(ub_swWIntake_mgd(i), temp_storage + ub_infil_mgd(i)); % assume as much SW as possible goes into the spreading basin
    for j = 1:numWRFs        
        if j == 1
            rwProduction0_cell{j}(i) = min(ub_rw1_mgd(i), (storageCapacity_mg - temp_storage) + ub_infil_mgd(i) - swInput0(i)); % any additional unused capacity can be used for RW
        else
            rwProduction0_cell{j}(i) = min(ub_rw2_mgd(i), (storageCapacity_mg - temp_storage) + ub_infil_mgd(i) - swInput0(i) - rwProduction0_cell{j}(i));
        end
    end
    
    infiltration0(i) = max(0, min(swInput0(i) + rwProduction1_0(i) + rwProduction2_0(i) + temp_storage - ET_mgd(i), ub_infil_mgd(i)));
    storage0(i) = max(0, temp_storage + swInput0(i) + rwProduction1_0(i) + rwProduction2_0(i) - infiltration0(i) - ET_mgd(i));
    
end

rwProduction0_array = [rwProduction0_cell{1}, rwProduction0_cell{2}];

x0 = lb;
disp('using lower bound as initial solution guess')