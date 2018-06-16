% Copyright Jonathan L. Bradshaw 2018

function [f, df] = funit_lcc_daily_nonlinearCost_twoWRFserial(z, rwCapacity_mgd,a,b,c,d,onlineFactor, om_fixed_treat_fraction,assessmentPeriod_yr,conveyCost_fixed,H_ft,hfs_perQ_mgd,energyCost_perMGDperftperday, alpha_param, rw_indices, infiltration_indices, discountFactor_average)

% To include the pipeline cost, this updated function includes the lifecycle pipeline
% construction & O&M cost exlcuding energy (pipeCost_constrOM), the net elevation head loss (H), and
% the frictional head loss without the flow term (hfs_perQ_afy = hf /
% Q^1.85), the life cycle cost factor for energy (energyLCCFactor), the
% cost of energy in $/kWh (energyCost_USDperkWH)

% Purpose: anonymous function that allows for passing extra parameters into
% the optimization command.
% http://www.mathworks.com/help/optim/ug/passing-extra-parameters.html
% Compute the unit life cycle cost (i.e., life cycle cost per
% volume nominal water recharged) and Jacobian (i.e., gradient vector)

% Output: Value of unit life cycle cost (f) and gradient vector (df) given
% the inputs.

numWRFs = length(rwCapacity_mgd);

alpha_pattern = zeros(1,numWRFs);
alpha1 = alpha_pattern;
alpha2 = alpha_pattern;
alpha3 = alpha_pattern;
alpha4 = alpha_pattern;

if numWRFs > 2
    disp('Need to use another method for alpha3')
end

x_mgd = sum(rwCapacity_mgd); % For the production costs, capacity is based on the sum of the two plants

for i = 1:numWRFs
    
    capitalCost_treat = b*x_mgd^a;
    om_treat_fixedCost = (discountFactor_average*assessmentPeriod_yr)*om_fixed_treat_fraction*d*x_mgd^c;
    
    % The alpha 4 term represents the unit variable cost (e.g., variable
    % treatment costs = z(13)*OM_var_unit)
    om_treat_varCost = discountFactor_average*(1-om_fixed_treat_fraction)*d*x_mgd^c/(x_mgd*onlineFactor)/units_check(1,'year2days'); % The variable cost is in
    
    alpha2(i) = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd(i); % constant that coefficent on (daily flow)^2.85 term
    
    if i == 1
        alpha1(i) = capitalCost_treat + om_treat_fixedCost + conveyCost_fixed(i); % Alpha 1 (fixed costs) is constant for fixed capacity
        alpha3(i) = sum((z(rw_indices{1}) + z(rw_indices{2})).^2.85); % variable
        alpha4(i) = om_treat_varCost + discountFactor_average*energyCost_perMGDperftperday*H_ft(i); % constant that is coefficient to linear daily flow term
    else
        alpha1(i) = conveyCost_fixed(i);
        alpha3(i) = sum(z(rw_indices{i}).^2.85); % variable
        alpha4(i) = discountFactor_average*energyCost_perMGDperftperday*H_ft(i); % The energy cost to overcome the elevation head difference. This needs to be multiplied by z(26)
    end
    
end

total_rw_mg = sum(z(rw_indices{1}) + z(rw_indices{2})); % variable
total_infil_mg = sum(z(infiltration_indices)); % variable

f = sum(alpha1) + dot(alpha2,alpha3) + alpha4(1)*total_rw_mg + alpha4(2)*sum(z(rw_indices{2}))- alpha_param*discountFactor_average*total_infil_mg;

%     need to use chain rule here for derivative at terminal treatment
%     plant
df_dz_treat1 = @(i) alpha2(1)*2.85*(z(rw_indices{1}(i)) + z(rw_indices{2}(i)))^1.85 + alpha4(1); % derivative of the function with respect to production at termainal plant
df_dz_treat2 = @(i) df_dz_treat1(i) + alpha2(2)*2.85*z(rw_indices{2}(i))^1.85 + alpha4(2); % derivative of the function with respect to production at upstream plant
df_dz_infil = (-alpha_param*discountFactor_average);

df = zeros(1,length(z));

for i = 1:length(infiltration_indices)
    df(rw_indices{1}(i)) = df_dz_treat1(i);
    df(rw_indices{2}(i)) = df_dz_treat2(i);
    df(infiltration_indices(i)) = df_dz_infil;
end

