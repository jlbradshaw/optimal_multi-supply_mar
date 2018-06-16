% Copyright Jonathan L. Bradshaw 2018

function [minULCC_array, xstar_array_postProcess, energyCost_nominal, ET_actual_mgd, varCost_1mg_array, problemCaps] = opt_paper2(elev_profile_m, rw_capacityVector_mgd, a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat_fraction, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailable_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_potential_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, min_rwCapacity_mgd, max_rwCapacity_mgd, assessmentPeriod_yr, discountFactor_average)

%% initialize and define variables

numCapacityPoints = length(rw_capacityVector_mgd);

minULCC_array = zeros(numCapacityPoints,1); % solution: array of minimum ULCCs
xstar_array_raw = zeros(numCapacityPoints,numberDecisionVariables); % solution array
xstar_array_postProcess = zeros(numCapacityPoints,numberDecisionVariables); % solution after post-processing
energyCost_nominal = minULCC_array; % array for total conveyance energy cost
varCost_1mg_array = minULCC_array;
ET_actual_mgd = zeros(numCapacityPoints, number_timesteps);
problemCaps = cell(1,numCapacityPoints);

% Get the LCC factors
[~, ~, lccFactor_constr_convey_pump, lccFactor_constr_convey_pipe, ~, ~] = getLccFactor(assessmentPeriod_yr);


%% Constraints

% Import the pipe options
po = pipe_op;
po = import_pipe_types(po, 'pipetemplate.csv');

z_solutions = NaN(numCapacityPoints,numberDecisionVariables,length(po.diameters));
f_solutions = NaN(numCapacityPoints,length(po.diameters));

energyCost_solutions = f_solutions; % create a matrix with same structure as f_solutions that will store the energy cost of the system
varCost_1mg_solutions = f_solutions;
output_solutions = cell(numCapacityPoints,1);

elev_profile_ft = elev_profile_m*units_check(1,'m2ft');


for i = 1:numCapacityPoints
        
    rwCapacity_mgd_point = rw_capacityVector_mgd(i);
    
    rwCapacity_afy = rwCapacity_mgd_point*units_check(1,'mgd2afy');
    [pipeLength_ft, pipeFeasible_index, pump_reqs, ~] = pipeline_routine(elev_profile_m, rwCapacity_afy, po);
    
    H_ft = elev_profile_ft(end,2) - elev_profile_ft(1,2); % overall elevation head difference between start and end of pipeline
    
    phi = 4.73; % assumes US units
    
    hfs_perQ_mgd = 1.1*pipeLength_ft*phi/units_check(1,'cfs2mgd')^1.85./po.Cs(pipeFeasible_index).^1.85./po.diameters(pipeFeasible_index).^4.87; % Compute the head loss in the pipe excluding the flow term. The leading 1.1 coefficient represents the minor fractional head loss of 10%.    
    
    % For each feasible pipe size, find the minimum unit life cycle cost configuration
    pipeFeasible_indexNum = find(pipeFeasible_index);
    output_solutions{i} = cell(length(pipeFeasible_indexNum),1);

    [A_matrix, b_vector, Aeq, beq, lb, ub, rw_indices, ~, infiltration_indices, storage_indices, x0] = getConstraints_papertwo_singleWRF(rwCapacity_mgd_point, number_timesteps, wasteFractionMFRO, ET_potential_mgd, wwAvailable_mgd, swAvailable_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, storageCapacity_mg, infilCapacity_mgd);
    nonlcon = [];
    
    for j = 1:length(pipeFeasible_indexNum) % iternate over all the feasible pipe sizes
        
        % Compute the fixed costs for pipeline and pump stations--note:
        % this assumes pipe costs are already in 2015 dollars.
        pipeCost_cap = pipeLength_ft*(po.unit_cost_const(pipeFeasible_indexNum(j))*lccFactor_constr_convey_pipe);
        pipeCost_om = pipeLength_ft*(po.unit_cost_OM(pipeFeasible_indexNum(j))*discountFactor_average*assessmentPeriod_yr);
        pipeCost_fixed = pipeCost_cap + pipeCost_om; % Pipeline fixed costs
        
        pumpCostConstr = USD2011to2015(1)*3.12*10^(0.7583*log10(rwCapacity_mgd_point*units_check(1,'mgd2gpm')) + 3.1951); % Life cycled construction cost of pump
        pumpCostCap = lccFactor_constr_convey_pump*pumpCostConstr; % Life cycle capital cost of pump
        pumpCostOM = (assessmentPeriod_yr*discountFactor_average)*(USD2011to2015(2)*1e4 + 0.05*pumpCostConstr); % life cycle cost of non-energy pump O&M; Updated September 2017
        numberPumps = length(pump_reqs{j}(:,1));
        
        conveyCost_fixed = pipeCost_fixed + numberPumps*(pumpCostCap + pumpCostOM);
        
        alpha_param=0;                                                                    %Auxiliary parameter for linearization
        fstar=10;                                                                %Dummy first value
        opt_precision=1e-5;
        
        iter = 0
        
        %         [i,j]
        
        while abs(fstar)>opt_precision                        
            if iter > 25
                if abs(fstar) < opt_precision*1e1 % OK if it's within 10x of the precision
                    problemCaps{i} = [j, fstar];                    
                    break
                elseif sum(full(x0 ~= lb)) > 0 % if x0 is not already equal to the lower bound, set it to the lower bound
                    x0 = lb;
                    iter = 0;                    
                    disp('resetting optimization starting point to lower bound')               
                else
                    rwCapacity_mgd_point
                    j
                    disp('Optimization algorithm may be stuck')
%                     break             
                end
            end
            
            if iter > 30 && abs(fstar) < opt_precision*1e2                                
                problemCaps{i} = [problemCaps{i}; [pipeFeasible_indexNum(j), fstar]];
                break                
            end
            
            if iter > 35 && abs(fstar) < opt_precision*1e3                                
                problemCaps{i} = [problemCaps{i}; [pipeFeasible_indexNum(j), fstar]];
                break                
            end
            
            if iter > 40 && abs(fstar) < opt_precision*5e3                                
                problemCaps{i} = [problemCaps{i}; [pipeFeasible_indexNum(j), fstar]];
                break                
            end
            
            f_df = @(z) funit_lcc_daily_nonlinearCost(z, rwCapacity_mgd_point,a,b,c,d,om_fixed_treat_fraction,assessmentPeriod_yr,conveyCost_fixed,H_ft,hfs_perQ_mgd(j),energyCost_perMGDperftperday,onlineFactor, alpha_param, rw_indices, infiltration_indices, discountFactor_average);
            
            h_anon_daily = @(z,lambda) hessianfcn_daily(z,lambda, hfs_perQ_mgd(j), energyCost_perMGDperftperday, rw_indices, discountFactor_average);
            
            options = optimoptions('fmincon','Algorithm','interior-point', 'GradObj','on','HessFcn',h_anon_daily,'Display','off');
            
            [xstar,fstar] = fmincon(f_df,x0,A_matrix,b_vector,Aeq,beq,lb,ub,nonlcon,options); % fstar is the LCC/(nominal lifetime water volume)
            
            total_infil_mg = sum(xstar(infiltration_indices));
            total_cost = fstar+alpha_param*discountFactor_average*total_infil_mg;
            
            alpha_param = total_cost/(total_infil_mg*discountFactor_average);
            
            fstar
            iter = iter+1
            
            
        end
        
        z_solutions(i,:,pipeFeasible_indexNum(j)) = full(xstar);
        
        f_solutions(i,pipeFeasible_indexNum(j)) = alpha_param;
        
        alpha2_nominal = energyCost_perMGDperftperday*hfs_perQ_mgd(j); % Note: this is nominal cost, not PV cost
        alpha2 = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd(j); % Note: this is nominal cost, not PV cost
        alpha3 = sum(xstar(rw_indices).^2.85); % variable
        
        x_mgd = rwCapacity_mgd_point;
                
        om_treat_varCost = discountFactor_average*(1-om_fixed_treat_fraction)*d*x_mgd^c/(x_mgd*onlineFactor)/units_check(1,'year2days'); % The variable cost is in
        alpha4 = om_treat_varCost + discountFactor_average*energyCost_perMGDperftperday*H_ft;
        
        energyCost_convey_total = alpha3 * alpha2_nominal + energyCost_perMGDperftperday * H_ft * sum(xstar(rw_indices)); % compute the total lifetime energy cost
        if energyCost_convey_total < 0
            disp('Warning: negative energy use for conveyance')
        end
        energyCost_solutions(i,pipeFeasible_indexNum(j)) = energyCost_convey_total;
        
        varCost_1mg_solutions(i, pipeFeasible_indexNum(j)) = 2.85*alpha2 + alpha4;                
    end
        
    [temp_optimum, temp_optimum_index] = min(f_solutions(i,:));
    
    minULCC_array(i) = temp_optimum;
    xstar_array_raw(i,:) = z_solutions(i,:,temp_optimum_index);
    energyCost_nominal(i) = energyCost_solutions(i, temp_optimum_index);
    varCost_1mg_array(i) = varCost_1mg_solutions(i, temp_optimum_index);

    [xstar, ET_out_mgd] = postProcess(xstar_array_raw(i,:));
    
    xstar_array_postProcess(i,:) = xstar;
            
    ET_actual_mgd(i,:) = ET_out_mgd;
                  
end


    function [xstar_out, ET_out_mgd] = postProcess(xstar_in)
        
        % Purpose: correct storage and infiltration and ET, which may not be correct if
        % there is no driver to have infiltration occur as soon as possible
        
        %% Correct infiltration (i.e., use storage to exhause infiltration capacity)
        infilCapacity_mgd_array = ub(infiltration_indices);
        
        xstar_out = xstar_in;
        
        for k = 1:number_timesteps
            % Get the storage and infil values
            temp_storage1 = xstar_out(storage_indices(k));
            temp_infil = xstar_out(infiltration_indices(k));                                    
            
            temp_unused_infil = infilCapacity_mgd_array(k) - temp_infil; % the additional unused percolation capacity available
            
            if temp_unused_infil > 0 % if there could be more infiltration
                if temp_storage1 > temp_unused_infil % if there is more storage than could be infiltrated
                    temp_infil_delta = temp_unused_infil; % then increase infiltration by the remaining unused capacity
                else % else all the stored water can be infiltrated
                    temp_infil_delta = temp_storage1;
                end
            else
                temp_infil_delta = 0;
            end
            
            % Update infiltration and storage with these changes
            xstar_out(infiltration_indices(k)) = xstar_out(infiltration_indices(k)) + temp_infil_delta;
            xstar_out(storage_indices(k)) = xstar_out(storage_indices(k)) - temp_infil_delta;
            
            % Also modify xstar_out(infiltration_indices(k+1)) because of
            % the affect on storage(i)
            if k < number_timesteps
                xstar_out(infiltration_indices(k+1)) = xstar_out(infiltration_indices(k+1)) - temp_infil_delta; % note: do not need to use a max(0,...) here because we know that this is greater than 0 because infil(t) >= storage(t), and storage(t) above is modified such that it must be greater than or equal to 0.
            end
        end
        
        %% Correct storage (i.e., remove ET from storage)        
        ET_out_mgd = zeros(1,length(ET_potential_mgd));
        
        for k = 1:number_timesteps
            
            % Get the current storage and potential ET values
            temp_storage1 = xstar_out(storage_indices(k));
            temp_ET = ET_potential_mgd(k);
            
            if temp_storage1 > 0
                delta_ET = min(temp_ET, temp_storage1); % compute the actual ET based on the volume of water in storage. i.e., ET is equal to storage up to potential ET value)
                ET_out_mgd(k) = delta_ET; % assign the actual ET
                
                % Overwrite the affected storage(i) and infiltration(i+1) values
                % with the ET affect
                xstar_out(storage_indices(k)) = xstar_out(storage_indices(k)) - delta_ET;
                
                if k < number_timesteps
                    xstar_out(infiltration_indices(k+1)) = xstar_out(infiltration_indices(k+1)) - delta_ET; % I think a max(0, ...) condition is not required because infiltration >= storage, which must be >= 0 based on the above commands.
                end
                
            end           
            
        end        
              
    end




end

