% Copyright Jonathan L. Bradshaw 2018

function [f, df] = funit_lcc_daily_nonlinearCost(z, rwCapacity_mgd,a,b,c,d,om_fixed_treat_fraction,assessmentPeriod_yr,conveyCost_fixed,H_ft,hfs_perQ_mgd,energyCost_perMGDperftperday,onlineFactor, alpha_param, rw_indices, infiltration_indices, discountFactor_average)

% Purpose: anonymous function that allows for passing extra parameters into
% the optimization command.
% http://www.mathworks.com/help/optim/ug/passing-extra-parameters.html
% Compute the unit life cycle cost (i.e., life cycle cost per
% volume nominal water recharged) and Jacobian (i.e., gradient vector)

% Input: vector z (decision variables); rwCapacity_MGD is the plant capacity in MGD; a and b
% are the cost curve parameters for lifecycle construction cost; c and d
% are cost curve parameters for lifecycle O&M costs; om_fixted_treat is the fraction of O&M costs that are fixed; lifetime_yr is the
% lifetime of the project. Additional conveyance related inputs: pipeCost_constrOM is the
% lifecycle construction and O&M (exlcuding energy) cost for a pipeline;
% H_ft is the net head difference in ft between the starting and ending
% points on the pipeline; hfs_perQ is the equation for the friction and
% minor head losses in a pipeline omitting the Q terms.

% Output: Value of unit life cycle cost (f) and gradient vector (df) given
% the inputs.
   
x_mgd = rwCapacity_mgd;
capitalCost_treat = b*x_mgd^a;
om_treat_fixedCost = (discountFactor_average*assessmentPeriod_yr)*om_fixed_treat_fraction*d*x_mgd^c;
alpha1 = capitalCost_treat + om_treat_fixedCost + conveyCost_fixed; % Alpha (fixed costs) 1 is constant for fixed capacity 
alpha2 = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd; % constant
alpha3 = sum(z(rw_indices).^2.85); % variable
total_rw_mg = sum(z(rw_indices)); % variable
total_infil_mg = sum(z(infiltration_indices)); % variable

% The alpha 4 term represents the unit variable cost (e.g., variable
% treatment costs = z(13)*OM_var_unit)
om_treat_varCost = discountFactor_average*(1-om_fixed_treat_fraction)*d*x_mgd^c/(x_mgd*onlineFactor)/units_check(1,'year2days'); % The variable cost is in  
alpha4 = om_treat_varCost + discountFactor_average*energyCost_perMGDperftperday*H_ft; 

f = alpha1 + alpha2*alpha3 + alpha4*total_rw_mg - alpha_param*discountFactor_average*total_infil_mg;

df_dz_treat = @(zi) (alpha2*2.85*(zi)^1.85 + alpha4); % derivative of the function with respect to a specific month in afm

df_dz_infil = @(zi) (-alpha_param*discountFactor_average);

df = zeros(1,length(z));
for i = 1:length(rw_indices)
    df(rw_indices(i)) = df_dz_treat(z(rw_indices(i)));
end

for i = 1:length(infiltration_indices)
    df(infiltration_indices(i)) = df_dz_infil(z(infiltration_indices(i)));
end

end