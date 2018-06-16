% Copyright Jonathan L. Bradshaw 2018

function [capFactor,OMFactor] = computeLccFactor(discEff,lifetime,assessPeriod,fracComponent,con2capFactor,isLoan,loanRate,payback,discNom)

%% Input: 
% disc = effective discount rate (i.e., the net of nominal discount rate and escalation rate) (decimal)
% lifetime = vector with lifetime of components (years)
% assessPeriod = project assessment period (years)
% fracComponent = vector with fraction of capital cost for component (dimensionless)
% con2capFactor = contingency and implementation factor that converts
% construction to capital
% isLoan = boolean that indicates financing: 0 if it is pay-as-you-go (no loan), 1 if it is a loan
% loanRate = loan interest rate (decimal)
% payback = Payback period (years)
% discNom = nominal discount rate (decimal); note: only used for computing
% the present value of loan payments. 


numComponent = length(lifetime);

% Compute the annual payment assuming 1$ of capital

if isLoan==1
    A_perCap = loanRate*(1+loanRate)^payback/((1+loanRate)^payback-1);
          
    cap_PVik = fracComponent * A_perCap/ discNom *(1-1/(1+discNom)^payback);
    
else cap_PVik = fracComponent;    
end

rep_i_vector = zeros(assessPeriod*numComponent,1); % initalize vector to store the replacement costs
counter = 1;
for k = 1:numComponent % iterate over the components
    N = floor((assessPeriod-1)/lifetime(k));
    if N>0;
        for n = 1:floor((assessPeriod-1)/lifetime(k)) % iterate over the number of replacements
            rep_i_vector(counter) = cap_PVik(k)*(1/(1+discEff))^(n*lifetime(k));
            counter = counter + 1;            
        end
    end
end

rep_i = sum(rep_i_vector); % compute the total residual value factor per $1 of capital cost

% compute the salvage value per $1 capital
sal_ik = zeros(2,1);
for k = 1:numComponent
% salvage value is discounted based on the the end of the assessment
% period, not when it was constructed.
    sal_ik(k) = 1/con2capFactor * fracComponent(k) * (1/(1+discEff))^(assessPeriod)* (1-mod(assessPeriod,lifetime(k))/lifetime(k));
end

sal_i = sum(sal_ik);

capFactor = sum(cap_PVik)+rep_i-sal_i; % the capital multiplying factor includes the original capital costs plus replacement costs, minus salvage costs

OMFactor = sum((1/(1+discEff)).^(0:(assessPeriod-1))); % sum of the escalated and discounted values over the assessment period per $1 of O&M cost in the first year

sum(cap_PVik);