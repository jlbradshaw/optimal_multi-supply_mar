% Copyright Jonathan L. Bradshaw 2018

function Hout = hessianfcn_daily(z,lambda,hfs_perQ_mgd,energyCost_perMGDperftperday, rw_indices, discountFactor_average)

H = sparse(length(z), length(z));        
alpha2 = discountFactor_average*energyCost_perMGDperftperday*hfs_perQ_mgd; % constant

df2_dz_treat2 = @(zi) (alpha2*2.85*1.85*(zi)^0.85); 

for i = 1:length(rw_indices)
    H(rw_indices(i),rw_indices(i)) = df2_dz_treat2(z(rw_indices(i)));
end       
Hout = H;