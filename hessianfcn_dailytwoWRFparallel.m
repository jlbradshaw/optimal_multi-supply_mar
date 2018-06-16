% Copyright Jonathan L. Bradshaw 2018

function Hout = hessianfcn_dailytwoWRFparallel(z,lambda,rwCapacity_mgd,~,~,c,d,onlineFactor,om_fixed_treat_fraction,~,~,H_ft,hfs_perQ_mgd,energyCost_perMGDperftperday, rw_indices, infiltration_indices, discountFactor_average)

H = sparse(length(z), length(z));

numWRF = length(rwCapacity_mgd);

alpha_pattern = zeros(1,numWRF);
alpha2 = alpha_pattern;
om_treat_fixedCost = alpha_pattern;
x_mgd = alpha_pattern;
capitalCost_treat = alpha_pattern;
om_treat_varCost = alpha_pattern;

if numWRF > 2
    disp('Need to use another method for alpha3')
end

for i = 1:numWRF
    
    x_mgd(i) = rwCapacity_mgd(i);     
    
    om_treat_varCost(i) = discountFactor_average*(1-om_fixed_treat_fraction)*d*x_mgd(i)^c/(x_mgd(i)*onlineFactor)/units_check(1,'year2days'); % The variable cost is in       
    
    alpha2(i) = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd(i); % constant that coefficent on (daily flow)^2.85 term

    
end

df2_drw2 = @(wrf,t) (alpha2(wrf)*2.85*1.85*(z(rw_indices{wrf}(t)))^0.85);

for i = 1:length(infiltration_indices) % loop through the number of timesteps (equal to length of infiltration_indices)
    for j = 1:numWRF
        H(rw_indices{j}(i), rw_indices{j}(i)) = df2_drw2(j,i);        
    end
end
Hout = H;