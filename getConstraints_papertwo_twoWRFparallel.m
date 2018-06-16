% Copyright Jonathan L. Bradshaw 2018

function [A, b, Aeq, beq, lb, ub, rw1_indices, rw2_indices, sw_indices, infiltration_indices, storage_indices, x0] = getConstraints_papertwo_twoWRFparallel(rw1Capacity_mgd_point, rw2Capacity_mgd_point, number_timesteps, wasteFractionMFRO, ET_mgd, ww1Available_mgd, ww2Available_mgd, swAvailable_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, min_rw1Capacity_mgd, max_rw2Capacity_mgd, min_rw2Capacity_mgd, max_rw2Capacity_mgd, storageCapacity_mg, infilCapacity_mgd)

%% Indices
% Process capacities
ind_CFAT1=1;
ind_CFAT2=2;
number_processes = 1;                                                        % total number of processes represented

% Flows
number_flows = 2;
flow_cell = cell(1,number_flows);           % Create a 1-D array to store the flow indices


% loop through number of flows and create indices to assign to the
% different flows.
for i = 1:number_flows
    flow_cell{i} = (1:1:number_timesteps)+(i-1)*number_timesteps + number_processes;
end
rw1_indices = flow_cell{1};
rw2_indices = flow_cell{2};

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

ind_rw1_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,rw1_indices);                    %Creates linear index for q0
ind_rw2_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,rw2_indices);                    %Creates linear index for q0

ind_sw_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,sw_indices);
ind_infil_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,infiltration_indices);
ind_storage_A=sub2ind(size(A_timestepTemplate),1:number_timesteps,storage_indices);


%% Constraints
% Condition 1: incoming flow must not exceed capacity
A1_1 = A_timestepTemplate;
A1_1(1:number_timesteps,ind_CFAT1)=-1;                                             %Cbac coefficient for first condition
A1_1(ind_rw1_A)=1;                                                                  %q0 coefficient for first condition
b1_1=zeros(number_timesteps,1);

A1_2 = A_timestepTemplate;
A1_2(1:number_timesteps,ind_CFAT2)=-1;                                             %Cbac coefficient for first condition
A1_2(ind_rw2_A)=1;                                                                  %q0 coefficient for first condition
b1_2=zeros(number_timesteps,1);

A1 = [A1_1; A1_2];
b1 = [b1_1; b1_2];

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

%% Bounds
lb_rw1_mgd = zeros(1,number_timesteps);
lb_rw2_mgd = zeros(1,number_timesteps);
lb_sw_mgd = zeros(1,length(sw_indices));
lb_infil_mgd = zeros(1,length(infiltration_indices));
lb_storage_mg = zeros(1,length(storage_indices));

ub_storage_mg = storageCapacity_mg*ones(1,number_timesteps);
ub_infil_mgd = infilCapacity_mgd*ones(1,number_timesteps);
ub_swWIntake_mgd = min(swAvailable_mgd',intakeLimits_mgd);

A9_1=zeros(1,number_decisionVars);
A9_1(rw1_indices)=-1;
b9_1=-min_prodRW_fractionofCap*rw1Capacity_mgd_point*number_timesteps;

A9_2=zeros(1,number_decisionVars);
A9_2(rw2_indices)=-1;
b9_2=-min_prodRW_fractionofCap*rwCapacity_mgd_point2*number_timesteps;

A9 = [A9_1; A9_2];
b9 = [b9_1; b9_2];

% Compile all the constraints

A = sparse([A1; A9]);
b = sparse([b1; b9]);

Aeq = sparse(Aeq1);
beq = sparse(beq1);

% Specify additional bounds
lb_fat1Capacity_mgd = rw1Capacity_mgd_point;
ub_fat1Capacity_mgd = rw1Capacity_mgd_point;
ub_rw1_mgd = min(ww1Available_mgd, ub_fat1Capacity_mgd);

lb_fat2Capacity_mgd = rw2Capacity_mgd_point;
ub_fat2Capacity_mgd = rw2Capacity_mgd_point;
ub_rw2_mgd = min(ww2Available_mgd, ub_fat2Capacity_mgd);

if is_uncap_fixed % if the unused capacity is supposed to be fixed (e.g., by month)
    
    % replace the upper bound of the rw input with the alternative unused capacity
    ub_rw1_mgd = min(uncap_daily_mgd, ww1Available_mgd);
    ub_rw2_mgd = min(uncap_daily_mgd, ww2Available_mgd);
    
    if ~is_minRWCapacity_fixed
        disp('need to update this constraint, especially Aeq1'
    end
        
end

% Compile the lower and upper bounds
lb = sparse([lb_fat1Capacity_mgd, lb_fat2Capacity_mgd, lb_rw1_mgd, lb_rw2_mgd, lb_sw_mgd, lb_infil_mgd, lb_storage_mg]);
ub = sparse([ub_fat1Capacity_mgd, ub_fat2Capacity_mgd, ub_rw1_mgd, ub_rw2_mgd, ub_swWIntake_mgd, ub_infil_mgd, ub_storage_mg]);

%% Compute an initial guess for interior point method

rwCapacity1_0 = ub_fat1Capacity_mgd; % not an important guess because capacity is fixed
rwCapacity2_0 = ub_fat2Capacity_mgd;
swInput0 = zeros(1,number_timesteps);
storage0 = zeros(1,number_timesteps);
rwProduction1_0 = zeros(1,number_timesteps);
rwProduction2_0 = zeros(1,number_timesteps);
infiltration0 = zeros(1,number_timesteps);

initial_storage = 0;

for i = 1:number_timesteps
    if i == 1
        temp_storage = initial_storage;
    else
        temp_storage = storage0(i-1);
    end
    
    swInput0(i) = min(ub_swWIntake_mgd(i), temp_storage + ub_infil_mgd(i)); % assume as much SW as possible goes into the spreading basin
    rwProduction1_0(i) = min(ub_rw1_mgd(i), temp_storage + ub_infil_mgd(i) - swInput0(i)); % any additional unused capacity can be used for RW
    rwProduction2_0(i) = min(ub_rw2_mgd(i), temp_storage + ub_infil_mgd(i) - swInput0(i) - rwProduction1_0(i)); % any additional unused capacity can be used for RW
    
    infiltration0(i) = max(0, min(swInput0(i) + rwProduction1_0(i) + rwProduction2_0(i) + storage0(i)- ET_mgd(i), ub_infil_mgd(i)));
    storage0(i) = temp_storage + swInput0(i) + rwProduction1_0(i) + rwProduction2_0(i) - infiltration0(i) + - ET_mgd(i);    
end

x0 = [rwCapacity1_0, rwCapacity2_0, rwProduction1_0, rwProduction2_0, swInput0, infiltration0, storage0];