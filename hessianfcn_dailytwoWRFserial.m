% Copyright Jonathan L. Bradshaw 2018

function Hout = hessianfcn_dailytwoWRFserial(z,lambda,rwCapacity_mgd,~,~,c,d,onlineFactor,om_fixed_treat_fraction,~,~,H_ft,hfs_perQ_mgd,energyCost_perMGDperftperday, rw_indices, infiltration_indices, discountFactor_average)

H = sparse(length(z), length(z));

numWRFs = length(rwCapacity_mgd);

alpha_pattern = zeros(1,numWRFs);
alpha2 = alpha_pattern;
alpha4 = alpha_pattern;

x_mgd = sum(rwCapacity_mgd); % For the production costs, capacity is based on the sum of the two plants

om_treat_varCost = discountFactor_average*(1-om_fixed_treat_fraction)*d*x_mgd^c/(x_mgd*onlineFactor)/units_check(1,'year2days'); % The variable cost is in

if numWRFs > 2
    disp('Need to use another method for alpha3')
end

for i = 1:numWRFs
    
    alpha2(i) = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd(i); % constant
    
    if i == 1
        alpha4(i) = om_treat_varCost + discountFactor_average*energyCost_perMGDperftperday*H_ft(i); % constant that is coefficient to linear daily flow term
    else
        alpha4(i) = discountFactor_average*energyCost_perMGDperftperday*H_ft(i); % The energy cost to overcome the elevation head difference. This needs to be multiplied by z(26)
    end
    
end

df2_drw1 = @(i) (alpha2(1)*2.85*1.85*(z(rw_indices{1}(i)) + z(rw_indices{2}(i)))^0.85);
df2_drw2 = @(i) (alpha2(1)*2.85*1.85*(z(rw_indices{1}(i)) + z(rw_indices{2}(i)))^0.85 + alpha2(2)*2.85*1.85*(z(rw_indices{2}(i)))^0.85); % for second derivative wrt upstream flow term only

for i = 1:length(infiltration_indices) % loop through the number of timesteps (equal to length of infiltration_indices)
    H(rw_indices{1}(i), rw_indices{1}(i)) = df2_drw1(i);
    H(rw_indices{1}(i), rw_indices{2}(i)) = df2_drw1(i);
    H(rw_indices{2}(i), rw_indices{1}(i)) = df2_drw1(i);
    H(rw_indices{2}(i), rw_indices{2}(i)) = df2_drw2(i);    
end

Hout = H;