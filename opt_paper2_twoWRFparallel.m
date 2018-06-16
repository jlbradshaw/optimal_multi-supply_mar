% Copyright Jonathan L. Bradshaw 2018

function [minULCC_array, xstar_array_postProcess, energyCost_solutionsNominal, ET_actual_mgd, problemCaps] = opt_paper2_twoWRFparallel(elev_profile_m, rw_capacityVector_mgd, a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailable_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_potential_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, assessmentPeriod_yr, discountFactor_average, filename)

filename_here = strcat(filename,'_',num2str(rw_capacityVector_mgd{2}),'MGDcapacity');

% Number of WRFs considered:
numWRFs = 2;

%% initialize and define variables

elev_profile_ft = elev_profile_m; % initialize cell array

elev_profile_ft = elev_profile_m; % initialize cell array

% if want to develop a vector of potential capacities:
for i = 1:numWRFs
    elev_profile_ft{i} = elev_profile_m{i}*units_check(1,'m2ft');
end

% Get the LCC factors
[~, ~, lccFactor_constr_convey_pump, lccFactor_constr_convey_pipe, lccFactor_om_convey, ~] = getLccFactor(assessmentPeriod_yr);

% tic

numCapacityPoints = zeros(1,2);
numCapacityPoints(1) = length(rw_capacityVector_mgd{1});
numCapacityPoints(2) = length(rw_capacityVector_mgd{2});

minULCC_array = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution: array of minimum ULCCs
xstar_array_raw = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution array
xstar_array_postProcess = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution after post-processing
energyCost_solutionsNominal = cell(numCapacityPoints(1),numCapacityPoints(2));
ET_actual_mgd = cell(numCapacityPoints(1),numCapacityPoints(2));
problemCaps = cell(numCapacityPoints(1),numCapacityPoints(2));

% Import the pipe options
po = pipe_op;
po = import_pipe_types(po, 'pipetemplate.csv');


%% instances when there is a pipeline:

if length(elev_profile_m{i}) > 1
    
    % For a single pipeline, it is easy enough to keep track of the
    % pipe sizes, but for multiple pipelines, it becomes more work than
    % I think is useful; consequently, optimum solutions for each pipe
    % size are not produced here.
    
    tic        
    
    for i = 1:numCapacityPoints(1) % use parfor line above for parrallel processing        
        for j = 1:numCapacityPoints(2)
            
            [f_out, x_out_raw, x_out_pp, energyCost_out, ET_out_mgd] = opt_subroutine;
            minULCC_array{i,j} = f_out;
            xstar_array_raw{i,j} = x_out_raw;
            xstar_array_postProcess{i,j} = x_out_pp;
            energyCost_solutionsNominal{i,j} = energyCost_out;
            ET_actual_mgd{i,j} = ET_out_mgd;
            
            if mod(j,10)==0
                [i,j]
                toc
            end            
            save(filename_here)
        end
    end        
    
    toc
    
end

    function [f_out, x_out_raw, x_out_pp, energyCostNominal_out, ET_actual_mgd] = opt_subroutine
                
        cellArrayTemplate = cell(numWRFs,1);
        pipeFeasible_index = cellArrayTemplate;
        pump_reqs = cellArrayTemplate;
        prv_reqs = cellArrayTemplate;
        pipeFeasible_indexNum = cellArrayTemplate;
        hfs_perQ_mgd = cellArrayTemplate;
        
        matTemplate = zeros(numWRFs,1); % matrix template
        H_ft = matTemplate;
        conveyCost_fixed = matTemplate;
        pipeLength_ft = matTemplate;
        pipeCost_fixed = matTemplate;
        pumpCostConstr = matTemplate;
        pumpCostOM = matTemplate;
        numberPumps = matTemplate;
        
        rwCapacity_mgd = matTemplate;
        rwCapacity_mgd(1) = rw_capacityVector_mgd{1}(i);
        rwCapacity_mgd(2) = rw_capacityVector_mgd{2}(j);
        rwCapacity_afy = rwCapacity_mgd*units_check(1,'mgd2afy');
        
        
        for m = 1:numWRFs
            
            [pipeLength_ft(m), pipeFeasible_index{m}, pump_reqs{m}, prv_reqs{m}] = pipeline_routine(elev_profile_m{m}, rwCapacity_afy(m), po);           
            
            H_ft(m) = elev_profile_ft{m}(end,2) - elev_profile_ft{m}(1,2); % overall elevation head difference between start and end of pipeline
            
            phi = 4.73; % assumes US units
            
            hfs_perQ_mgd{m} = 1.1*pipeLength_ft(m)*phi/units_check(1,'cfs2mgd')^1.85./po.Cs(pipeFeasible_index{m}).^1.85./po.diameters(pipeFeasible_index{m}).^4.87; % Compute the head loss in the pipe excluding the flow term. The leading 1.1 coefficient represents the minor fractional head loss of 10%.
            
            % For each feasible pipe size, find index numbers
            pipeFeasible_indexNum{m} = find(pipeFeasible_index{m});
            
        end                
        
        fstarMatrix = NaN*ones(length(pipeFeasible_indexNum{1}),length(pipeFeasible_indexNum{2}));
        xstarMatrix = NaN*ones(length(pipeFeasible_indexNum{1}),length(pipeFeasible_indexNum{2}),numberDecisionVariables);
        energyCostNominal_array = fstarMatrix; % initialize matrix with same structure as fstarMatrix to hold energy cost
        
        nonlcon = [];
        
        for k = 2:2 % based on the DCT alone runs, this is the optimal pipe size when for DCT capacity = 27 MGD            
            disp('DCT pipe size is fixed')                       
            
            temp_hf1 = hfs_perQ_mgd{1}(k);
            temp_hfs2 = hfs_perQ_mgd{2};                                               
            
            for m = 1:length(pipeFeasible_indexNum{2}) % iterate over all the feasible pipe sizes
                
                for n = 1:numWRFs
                    if n == 1
                        temp_index = k; % temp_index refers to the pipe size selection
                    else
                        temp_index = m;
                    end
                    
                    pipeCost_fixed(n) = pipeLength_ft(n)*(po.unit_cost_const(pipeFeasible_indexNum{n}(temp_index))*lccFactor_constr_convey_pipe + po.unit_cost_OM(pipeFeasible_indexNum{n}(temp_index))*lccFactor_om_convey); % Pipeline fixed costs                    
                    pumpCostConstr(n) = USD2011to2015(1)*lccFactor_constr_convey_pump*3.12*10^(0.7583*log10(units_check(rwCapacity_mgd(n),'mgd2gpm')) + 3.1951); % Life cycled construction cost of pump
                    pumpCostOM(n) = USD2011to2015(2)*lccFactor_om_convey*(1e4 + 0.05*pumpCostConstr(n)); % life cycle cost of non-energy pump O&M
                    numberPumps(n) = length(pump_reqs{n}{temp_index}(:,1));
                    
                    conveyCost_fixed(n) = pipeCost_fixed(n) + numberPumps(n)*(pumpCostConstr(n) + pumpCostOM(n));
                end
                
                hfs_vector = [temp_hf1, temp_hfs2(m)];
                
                [A_matrix, b_vector, Aeq, beq, lb, ub, rw_indices, ~, infiltration_indices, storage_indices, x0] = getConstraints_papertwo_twoWRFseries(rwCapacity_mgd, number_timesteps, wasteFractionMFRO, ET_potential_mgd, wwAvailable_mgd, swAvailable_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, storageCapacity_mg, infilCapacity_mgd);
                
                alpha_param=0;                                                                    %Auxiliary parameter for linearization
                fstar=10;                                                                %Dummy first value
                opt_precision=1e-3;
                
                iter = 0
                                
                
                while abs(fstar)>opt_precision
                    
                    if iter > 20
                        if abs(fstar) < opt_precision*5 % OK if it's within 10x of the precision
                            problemCaps{i,j} = [problemCaps{i,j}; pipeFeasible_indexNum{1}(k), pipeFeasible_indexNum{2}(m), fstar];
                            break
                        elseif sum(full(x0 ~= lb)) > 0 % if x0 is not already equal to the lower bound, set it to the lower bound
                            x0 = lb;
                            iter = 0;
                            disp('resetting optimization starting point to lower bound')
                        else
                            disp('Optimization algorithm may be stuck')
                            rwCapacity_mgd(1)
                            rwCapacity_mgd(2)
                            pipeFeasible_indexNum{1}(k)
                            pipeFeasible_indexNum{2}(m)
                            fstar
                        end
                    end
                                        
                    if iter > 25 && abs(fstar) < 0.01
                        problemCaps{i,j} = [problemCaps{i,j}; pipeFeasible_indexNum{1}(k), pipeFeasible_indexNum{2}(m), fstar];
                        disp('Quitting optimizaton: converged to solution within 1% after 25+ iterations. Check problemCaps for more information')
                        break
                    end
                    
                    if iter > 35 && abs(fstar) > 0.5                       
                        problemCaps{i,j} = [problemCaps{i,j}; pipeFeasible_indexNum{1}(k), pipeFeasible_indexNum{2}(m), fstar];
                        disp('Quitting optimizaton: could not converge to solution within 50% after 35+ iterations. Check problemCaps for more information')
                        break
                    end
                    
                    if iter > 40 && abs(fstar) > 0.05                       
                        problemCaps{i,j} = [problemCaps{i,j}; pipeFeasible_indexNum{1}(k), pipeFeasible_indexNum{2}(m), fstar];
                        disp('Quitting optimizaton: could not converge to solution within 5% after 40+ iterations. Check problemCaps for more information')
                        break
                    end
                    
                    if iter >= 50 && abs(fstar)                        
                        problemCaps{i,j} = [problemCaps{i,j}; pipeFeasible_indexNum{1}(k), pipeFeasible_indexNum{2}(m), fstar];
                        disp('Quitting optimizaton: Could not converge to solution within 50% after 50 iterations. Check problemCaps for more information')
                        break
                    end
                                        
                    f_df = @(z) funit_lcc_daily_nonlinearCost_twoWRFparallel(z, rwCapacity_mgd,a,b,c,d,onlineFactor,om_fixed_treat,assessmentPeriod_yr,conveyCost_fixed,H_ft,hfs_vector,energyCost_perMGDperftperday,alpha_param, rw_indices, infiltration_indices, discountFactor_average);
                    
                    hessianfcn = @(z,lambda) hessianfcn_dailytwoWRFparallel(z,lambda,rwCapacity_mgd,a,b,c,d,onlineFactor,om_fixed_treat,assessmentPeriod_yr,conveyCost_fixed,H_ft,hfs_vector,energyCost_perMGDperftperday, rw_indices, infiltration_indices, discountFactor_average);
                    

                    options = optimoptions('fmincon','Algorithm','interior-point', 'GradObj','on','HessFcn',hessianfcn,'Display','off');
                                        
                    [xstar,fstar] = fmincon(f_df,x0,A_matrix,b_vector,Aeq,beq,lb,ub,nonlcon,options); % fstar is the LCC/(nominal lifetime water volume)
                    
                    total_infil_mg = sum(xstar(infiltration_indices));
                    total_cost = fstar+alpha_param*discountFactor_average*total_infil_mg;
                    
                    alpha_param = total_cost/(total_infil_mg*discountFactor_average);
                    
                    fstar
                    iter = iter+1
                    
                    
                end
                
                disp('Exiting while loop in optimization function')
                                
                fstarMatrix(k,m) = alpha_param;
                
                xstarMatrix(k,m,:) = full(xstar);
                                
                
                %% Compute the conveyance energy use of this scenario (notation based on funit_lcc_monthly)
                alpha2 = energyCost_perMGDperftperday * hfs_vector;
                alpha3 = [sum(xstar(rw_indices{1}).^2.85), sum(xstar(rw_indices{2}).^2.85)];                
                
                energyCost_convey_total = dot(alpha3, alpha2) + energyCost_perMGDperftperday * (H_ft(1) * sum(xstar(rw_indices{1})) + H_ft(2) * sum(xstar(rw_indices{2}))); % compute the total lifetime energy cost. Note that the most downstream WRF, indexed as WRF 1 will have to pump the full load.
                if energyCost_convey_total < 0
                    disp('Warning: negative energy use for conveyance')
                end
                
                energyCostNominal_array(k,m) = energyCost_convey_total;
                                
            end
            
        end
        
        
        % Select the optimal f and z        
        f_out = min(min(fstarMatrix)); % Find the minimum value in the 3D matrix
        [row,col] = ind2sub(size(fstarMatrix),find(fstarMatrix==f_out)); % find the row, column, and page indices of the optimum solution within the fstarMatrix .
        x_out_raw = squeeze(xstarMatrix(row,col,:)); % produce the decision variable values associated with the optimum value of f.
        energyCostNominal_out = energyCostNominal_array(row,col);
                
        [x_out_pp, ET_actual_mgd] = postProcess(x_out_raw);        
               
        function [xstar_out, ET_out_mgd] = postProcess(xstar_in)
            
            % Purpose: correct storage and infiltration and ET, which may not be correct if
            % there is no driver to have infiltration occur as soon as possible                        
            
            
            %% Correct infiltration (i.e., use storage to exhaust infiltration capacity)
            infilCapacity_mgd_array = ub(infiltration_indices);
            
            xstar_out = xstar_in;
            
            for p = 1:number_timesteps
                % Get the storage and infil values
                temp_storage1 = xstar_out(storage_indices(p));
                temp_infil = xstar_out(infiltration_indices(p));
                
                temp_unused_infil = infilCapacity_mgd_array(p) - temp_infil; % the additional unused percolation capacity available
                
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
                xstar_out(infiltration_indices(p)) = xstar_out(infiltration_indices(p)) + temp_infil_delta;
                xstar_out(storage_indices(p)) = xstar_out(storage_indices(p)) - temp_infil_delta;
                
                % Also modify xstar_out(infiltration_indices(k+1)) because of
                % the affect on storage(i)
                if p < number_timesteps
                    xstar_out(infiltration_indices(p+1)) = xstar_out(infiltration_indices(p+1)) - temp_infil_delta; % note: do not need to use a max(0,...) here because we know that this is greater than 0 because infil(t) >= storage(t), and storage(t) above is modified such that it must be greater than or equal to 0.
                end
            end
            
            %% Then correct storage (i.e., remove ET from storage)
            ET_out_mgd = zeros(1,length(ET_potential_mgd));
            
            for p = 1:number_timesteps
                
                % Get the current storage and potential ET values
                temp_storage1 = xstar_out(storage_indices(p));
                temp_ET = ET_potential_mgd(p);
                
                if temp_storage1 > 0
                    delta_ET = min(temp_ET, temp_storage1); % compute the actual ET based on the volume of water in storage. i.e., ET is equal to storage up to potential ET value)
                    ET_out_mgd(p) = delta_ET; % assign the actual ET
                    
                    % Overwrite the affected storage(i) and infiltration(i+1) values
                    % with the ET effect
                    xstar_out(storage_indices(p)) = xstar_out(storage_indices(p)) - delta_ET;
                    
                    if p < number_timesteps
                        xstar_out(infiltration_indices(p+1)) = xstar_out(infiltration_indices(p+1)) - delta_ET; % I think a max(0, ...) condition is not required because infiltration >= storage, which must be >= 0 based on the above commands.
                    end
                    
                end                                                
                
            end
            
        end
        
        
    end


end



