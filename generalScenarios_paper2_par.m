% Copyright Jonathan L. Bradshaw 2018

function [] = generalScenarios_paper2_par(rw_capacityVector_mgd, assessmentPeriod_yr, scenario, is_uncap_fixed, ~, is_minRWCapacity_fixed, is_ETinConstraints, spreadingBasin_safetyFactor, filename)


% Purpose: Run the optimization process

%% Parameter definition

% Cost conversion using LA CCI
USD2011to2015_cci = 1.101;

% Cost conversion using LA CPI
USD2011to2015_cpi = 1.055; % This is from the California Consumer Price Index in LA-Anaheim area

% Create a vector that contains both cci and cpi cost conversion factors
USD2011to2015 = [USD2011to2015_cci, USD2011to2015_cpi];

energyCost_USDperkWh = 0.12;
onlineFactor = 0.92;


% Get the LCC factors
[lccFactor_constr_treat, lccFactor_om_treat, ~, ~, ~, ~] = getLccFactor(assessmentPeriod_yr);

discountFactor_average = lccFactor_om_treat/assessmentPeriod_yr; % factor to multiply annual costs by to account for average discounting

% Metric units
energyCost_USDperJ = energyCost_USDperkWh/units_check(1,'kwh2j');
density.water_kgperm3 = 1000;
pumpEfficiency = 0.75;
g.mpers2 = 9.81;
powerW_perm3sperm = density.water_kgperm3 * g.mpers2 / pumpEfficiency;
energyCost_perm3spermpers = energyCost_USDperJ*powerW_perm3sperm;
energyCost_perMGDperftperday = energyCost_perm3spermpers/units_check(1,'m3s2mgd')/units_check(1,'m2ft')/units_check(1,'s2d');

%% log-log cost model model: log_cost = a*log(capacity) + log(b); coefficients are for capacity in MGD

a_constr_treat = 0.8461;
lnb_constr_treat = 16.34;
b_constr_treat = exp(lnb_constr_treat);

% Operation costs
a_om_treat = 0.7970; % cost curve parameter
lnb_om_treat = 13.96; % cost curve parameter
b_om_treat = exp(lnb_om_treat); % cost curve parameter for MGD

om_fixed_treat = 0.42; % percentage of O&M costs that are fixed; computed in LA_cost_data spreadsheet
om_var_treat = 1 - om_fixed_treat; % percentage of O&M costs that are variable

a = a_constr_treat;
b = b_constr_treat*lccFactor_constr_treat;
c = a_om_treat;
d = b_om_treat; 

%% Read in pipeline information from ArcGIS data
route_filename = 'pipelineroutes.csv';

elev_prof_array_m = import_piperoutes(route_filename);

rwCapacity_DCT_mgd_permonth = 27*ones(1,12);
rwCapacity_LAG_mgd_permonth = [40*ones(1,6), 20*ones(1,6)];
rwCapacity_HTP_mgd_permonth = 160*ones(1,12);


%% Scenario general input parameters

%% Scenario general input parameters
% Time parameters
starttime = 1;
startyr = 1986;
t1 = datetime(startyr,10,1);
t2 = datetime(startyr+assessmentPeriod_yr, 10, 1);

% number_timesteps = floor(365*assessmentPeriod);
number_timesteps = caldays(between(t1,t2,'Days'));

numCapacityPoints = length(rw_capacityVector_mgd);
problemCaps = cell(1,numCapacityPoints);

% Recycled water parameters
wasteFractionMFRO = 0.2;                                                                %Percentage of waste from MF and RO
minPercentMFRO=0.05;%:0.05:0.95;                                                        %Minimum percentage of advanced treatment 5% to 95%

min_rwCapacity_mgd = 1;
max_rwCapacity_mgd = 100;

min_prodRW_fractionofCap = 1/100;                                                         % Average minimum water to be produced DAILY


%% Hansen parameters

% Stormater and spreading basin parameters
load('hansen_sw_availability_cfs_withDates.mat')
swAvailable_mgd = hansen_sw_availability_cfs_withDates(:,4)*units_check(1,'cfs2mgd');
swAvailable_mgd = swAvailable_mgd(starttime:number_timesteps + starttime - 1); %zeros(366,1);%                                 %Temporal, just pick a subset of the entire dataset


dateArray = hansen_sw_availability_cfs_withDates(:,1:3);
sw_available_monthArray = dateArray(:,2);
wateryear_months = [10:12, 1:9];

if is_uncap_fixed
    
    uncap_monthly_afm = [6904.258065	7246.8	5531.032258	3168.774194	3137.097345	0	0	3272.322581	4856.8	6516.193548	6472.645161	7330.8]; % monthly Unused capacity for Hansen (AF--normalized to 30 day month; computed in 'Spreding grounds rechrge.xls'
    uncap_monthly_mgd = units_check(1,'afm2mgd')*uncap_monthly_afm;
    
    uncap_daily_mgd = zeros(1,length(sw_available_monthArray));
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        uncap_daily_mgd(temp_index) = uncap_monthly_mgd(i);
    end
    
    uncap_daily_mgd = uncap_daily_mgd(starttime:number_timesteps + starttime - 1);
    
else
    % unimportant, but should define because we pass through the function
    uncap_monthly_mgd = NaN;
    uncap_daily_mgd = NaN;
    
end

% Evapotranspiration parameters
smallET = 0;
ET_mgd=smallET*ones(length(swAvailable_mgd),1);                                              %Negligible ET

potentialEvap_mgd = [0.640825124, 0.492392061, 0.387428423, 0.410989272, 0.489342875, 0.643852676, 0.792670632, 0.869216326, 1.00991129, 1.087478096, 1.048720028, 0.865205789];
for i=1:length(wateryear_months)
    temp_index = sw_available_monthArray == wateryear_months(i);
    ET_mgd(temp_index) = potentialEvap_mgd(i);
end

ET_mgd = ET_mgd(starttime:number_timesteps + starttime - 1);


storageCapacity_mg = 1400*units_check(1,'af2mg');                                                           %Max storage in Hansen: 1,400Acre-ft
storageCapacity_mg = storageCapacity_mg/spreadingBasin_safetyFactor;

infilCapacity_mgd = 300*units_check(1,'af2mg');                                                     %Percolation capacity rate in MGD
infilCapacity_mgd = infilCapacity_mgd/spreadingBasin_safetyFactor;


intakeCapacity_mgd = 600*units_check(1,'cfs2mgd');
intakeLimits_mgd = intakeCapacity_mgd*ones(1,number_timesteps);


%% Scenario: DCT to Hansen

if strcmp('DCT', scenario)
    
    numberDecisionVariables = number_timesteps*4;
    
    ev_m = elev_prof_array_m{1,2};
    
    rw_monthly = rwCapacity_DCT_mgd_permonth;
    
    maxRW = max(rwCapacity_DCT_mgd_permonth);
    
    wwAvailability_DCT_monthly_mgd = rw_monthly;
    wateryear_months = [10:12, 1:9];
    
    wwAvailability_DCT_daily_mgd = zeros(1,length(sw_available_monthArray));
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        wwAvailability_DCT_daily_mgd(temp_index) = wwAvailability_DCT_monthly_mgd(i);
    end
    
    wwAvailability_DCT_daily_mgd = wwAvailability_DCT_daily_mgd(starttime:number_timesteps + starttime - 1);
    
    wwAvailable_mgd = wwAvailability_DCT_daily_mgd;
        
    minULCC_DCTtoHansen_array = cell(1,numCapacityPoints);
    xstar_array_DCTtoHansen = cell(1,numCapacityPoints);
    energyCost_nominal_DCTtoHansen = cell(1,numCapacityPoints);
    ET_mgd_solutions_DCTtoHansen = cell(1,numCapacityPoints);
    varCost_1mg_array = cell(1,numCapacityPoints);        
    
    parfor i = 1:numCapacityPoints
%     for i = 1:numCapacityPoints
%         disp('not using parfor')        

        rw_capacityVector_mgd(i)
        [minULCC_DCTtoHansen, xstar, energyCost_nominal, ET_out, varCost_1mg, probCaps] = opt_paper2(ev_m, rw_capacityVector_mgd(i), a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailable_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, min_rwCapacity_mgd, maxRW, assessmentPeriod_yr, discountFactor_average);
        minULCC_DCTtoHansen_array{i} = minULCC_DCTtoHansen;
        xstar_array_DCTtoHansen{i} = xstar;
        energyCost_nominal_DCTtoHansen{i} = energyCost_nominal;
        ET_mgd_solutions_DCTtoHansen{i} = ET_out;
        varCost_1mg_array{i} = varCost_1mg;
        problemCaps{i} = probCaps{1};
        
    end    
    

%% Scenario: LAG to Hansen (with expansion)
elseif strcmp('LAG', scenario)
    
    numberDecisionVariables = number_timesteps*4;
    
    ev_m = elev_prof_array_m{9,2};
    
    rw_monthly = rwCapacity_LAG_mgd_permonth;
    
    maxRW = max(rwCapacity_LAG_mgd_permonth);
    
    wwAvailability_LAG_monthly_mgd = rw_monthly;
    wateryear_months = [10:12, 1:9];
    
    wwAvailability_LAG_daily_mgd = zeros(1,length(sw_available_monthArray));
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        wwAvailability_LAG_daily_mgd(temp_index) = wwAvailability_LAG_monthly_mgd(i);
    end
    
    wwAvailability_LAG_daily_mgd = wwAvailability_LAG_daily_mgd(starttime:number_timesteps + starttime - 1);
    
    wwAvailable_mgd = wwAvailability_LAG_daily_mgd;
    
    minULCC_LAGtoHansen_array = cell(1,numCapacityPoints);
    xstar_array_LAGtoHansen = cell(1,numCapacityPoints);
    energyCost_nominal_LAGtoHansen = cell(1,numCapacityPoints);
    ET_mgd_solutions_LAGtoHansen = cell(1,numCapacityPoints);
    varCost_1mg_array = cell(1,numCapacityPoints);
            
    parfor i = 1:numCapacityPoints        
        
        rw_capacityVector_mgd(i)
        
        [minULCC_LAGtoHansen, xstar, energyCost_nominal, ET_out, varCost_1mg, probCaps] = opt_paper2(ev_m, rw_capacityVector_mgd(i), a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailable_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, min_rwCapacity_mgd, maxRW, assessmentPeriod_yr, discountFactor_average);
        minULCC_LAGtoHansen_array{i} = minULCC_LAGtoHansen;
        xstar_array_LAGtoHansen{i} = xstar;
        energyCost_nominal_LAGtoHansen{i} = energyCost_nominal;
        ET_mgd_solutions_LAGtoHansen{i} = ET_out;
        varCost_1mg_array{i} = varCost_1mg;
        problemCaps{i} = probCaps{1};
        
    end

    
%% Scenario: Hyperion to Hansen

elseif strcmp('HTP', scenario)
    
    numberDecisionVariables = number_timesteps*4;
    
    ev_m = elev_prof_array_m{10,2};    
    
    rw_monthly = rwCapacity_HTP_mgd_permonth;
    
    maxRW = max(rwCapacity_HTP_mgd_permonth);
    
    wwAvailability_HTP_monthly_mgd = rw_monthly;
    wateryear_months = [10:12, 1:9];
    
    wwAvailability_HTP_daily_mgd = zeros(1,length(sw_available_monthArray));
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        wwAvailability_HTP_daily_mgd(temp_index) = wwAvailability_HTP_monthly_mgd(i);
    end
    
    wwAvailability_HTP_daily_mgd = wwAvailability_HTP_daily_mgd(starttime:number_timesteps + starttime - 1);
    
    wwAvailable_mgd = wwAvailability_HTP_daily_mgd;
        
    minULCC_HTPtoHansen_array = cell(1,numCapacityPoints);
    xstar_array_HTPtoHansen = cell(1,numCapacityPoints);
    energyCost_nominal_HTPtoHansen = cell(1,numCapacityPoints);
    ET_mgd_solutions_HTPtoHansen = cell(1,numCapacityPoints);
    varCost_1mg_array = cell(1,numCapacityPoints);
        
    
    counter = 1;
    numCap_remaining = numCapacityPoints;
    increment = min(10,numCapacityPoints);
    while numCap_remaining > 0
            
        for_vector = min(counter,numCapacityPoints):min(counter+increment,numCapacityPoints);
        parfor i = for_vector

            [minULCC_HTPtoHansen, xstar, energyCost_nominal, ET_out, varCost_1mg, probCaps] = opt_paper2(ev_m, rw_capacityVector_mgd(i), a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailable_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, min_prodRW_fractionofCap, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, min_rwCapacity_mgd, maxRW, assessmentPeriod_yr, discountFactor_average);
            minULCC_HTPtoHansen_array{i} = minULCC_HTPtoHansen;
            xstar_array_HTPtoHansen{i} = xstar;
            energyCost_nominal_HTPtoHansen{i} = energyCost_nominal;
            ET_mgd_solutions_HTPtoHansen{i} = ET_out;
            varCost_1mg_array{i} = varCost_1mg;
            problemCaps{i} = probCaps{1};
            
        end        
        save(filename)
        numCap_remaining = numCap_remaining - length(for_vector)       
        counter = counter + length(for_vector);       
    end    


%% Scenario 1c: LAG to DCT to Hansen (in series)

elseif strcmp('LAG-DCT_series', scenario)

    numberDecisionVariables = number_timesteps*5;
    
    % create the elevation profile cell array
    ev_m = cell(1,2);
    ev_m{1} = elev_prof_array_m{1,2}; % elevation profile of pipeline from DCT to Hansen
    ev_m{2} = elev_prof_array_m{3,2}; % elevation profile of pipeline from LAG to DCT
    
    % RW cell array
    rw_monthly = cell(1,2);
    rw_monthly{1} = rwCapacity_DCT_mgd_permonth;
    rw_monthly{2} = rwCapacity_LAG_mgd_permonth;
    
    % create arrays showing water availability per day
    wateryear_months = [10:12, 1:9];
    
    wwAvailability_daily_mgd = cell(1,2);
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        wwAvailability_daily_mgd{1}(temp_index) = rw_monthly{1}(i);
        wwAvailability_daily_mgd{2}(temp_index) = rw_monthly{2}(i);
    end
    
    % rewrite to only assign the days needed for the assessment period
    wwAvailability_daily_mgd{1} = wwAvailability_daily_mgd{1}(starttime:number_timesteps + starttime - 1);
    wwAvailability_daily_mgd{2} = wwAvailability_daily_mgd{2}(starttime:number_timesteps + starttime - 1);        
    
    numCapacityPoints = zeros(1,2);
    numCapacityPoints(1) = length(rw_capacityVector_mgd{1});
    numCapacityPoints(2) = length(rw_capacityVector_mgd{2});
    
    minULCC_DCTLAGtoHansen_array = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution: array of minimum ULCCs    
    xstar_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution after post-processing
    energyCost_nominal_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2));
    ET_mgd_solutions_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2));
    problemCaps = cell(numCapacityPoints(1),numCapacityPoints(2));
    
    if numCapacityPoints(1) > 1
        disp('This parfor method is not compatible for numCapacityPoints(1) > 1. May need to revert to method prior to 1/25/2018')
        keyboard
    end

parfor i = 1:numCapacityPoints(2)
    rw_capVector_temp = {rw_capacityVector_mgd{1}(1),rw_capacityVector_mgd{2}(i)};    
    [minULCC_DCTLAGtoHansen, xstar, energyCost_nominal, ET_out, probCaps] = opt_paper2_twoWRFseries(ev_m, rw_capVector_temp, a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailability_daily_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, assessmentPeriod_yr, discountFactor_average, filename);;
    minULCC_DCTLAGtoHansen_array{i} = minULCC_DCTLAGtoHansen{1};
    xstar_DCTLAGtoHansen{i} = xstar{1};
    energyCost_nominal_DCTLAGtoHansen{i} = energyCost_nominal{1};
    ET_mgd_solutions_DCTLAGtoHansen{i} = ET_out{1};
    problemCaps{i} = probCaps{1};
    
end
    save(filename)    
    plotScript_paper2(xstar_DCTLAGtoHansen, minULCC_DCTLAGtoHansen_array, ET_mgd_solutions_DCTLAGtoHansen, rw_indices, sw_indices, storage_indices, infiltration_indices, number_timesteps, rw_capacityVector_mgd, 2)

    
    elseif strcmp('LAG-DCT_parallel', scenario)    

    numberDecisionVariables = number_timesteps*5;
    
    % create the elevation profile cell array
    ev_m = cell(1,2);
    ev_m{1} = elev_prof_array_m{1,2}; % elevation profile of pipeline from DCT to Hansen
    ev_m{2} = elev_prof_array_m{9,2}; % elevation profile of pipeline from LAG to Hansen
    
    % RW cell array
    rw_monthly = cell(1,2);
    rw_monthly{1} = rwCapacity_DCT_mgd_permonth;
    rw_monthly{2} = rwCapacity_LAG_mgd_permonth;
    
    % create arrays showing water availability per day
    wateryear_months = [10:12, 1:9];
    
    wwAvailability_daily_mgd = cell(1,2);
    
    for i=1:length(wateryear_months)
        temp_index = sw_available_monthArray == wateryear_months(i);
        wwAvailability_daily_mgd{1}(temp_index) = rw_monthly{1}(i);
        wwAvailability_daily_mgd{2}(temp_index) = rw_monthly{2}(i);
    end
    
    % rewrite to only assign the days needed for the assessment period
    wwAvailability_daily_mgd{1} = wwAvailability_daily_mgd{1}(starttime:number_timesteps + starttime - 1);
    wwAvailability_daily_mgd{2} = wwAvailability_daily_mgd{2}(starttime:number_timesteps + starttime - 1);        
    
    numCapacityPoints = zeros(1,2);
    numCapacityPoints(1) = length(rw_capacityVector_mgd{1});
    numCapacityPoints(2) = length(rw_capacityVector_mgd{2});            
    
    minULCC_DCTLAGtoHansen_array = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution: array of minimum ULCCs    
    xstar_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2)); % solution after post-processing
    energyCost_nominal_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2));
    ET_mgd_solutions_DCTLAGtoHansen = cell(numCapacityPoints(1),numCapacityPoints(2));    
    problemCaps = cell(numCapacityPoints(1),numCapacityPoints(2));
    
    if numCapacityPoints(1) > 1
        disp('This parfor method is not compatible for numCapacityPoints(1) > 1. May need to revert to method prior to 1/25/2018')
        keyboard
    end

parfor i = 1:numCapacityPoints(2)
    
    rw_capVector_temp = {rw_capacityVector_mgd{1}(1),rw_capacityVector_mgd{2}(i)};

    [minULCC_DCTLAGtoHansen, xstar, energyCost_nominal, ET_out, probCaps] = opt_paper2_twoWRFparallel(ev_m, rw_capVector_temp, a, b, c, d, USD2011to2015, energyCost_perMGDperftperday, onlineFactor, om_fixed_treat, storageCapacity_mg, infilCapacity_mgd, swAvailable_mgd, wwAvailability_daily_mgd, intakeCapacity_mgd, numberDecisionVariables, number_timesteps, wasteFractionMFRO, ET_mgd, intakeLimits_mgd, uncap_monthly_mgd, uncap_daily_mgd, is_minRWCapacity_fixed, is_uncap_fixed, is_ETinConstraints, assessmentPeriod_yr, discountFactor_average, filename);;
    minULCC_DCTLAGtoHansen_array{i} = minULCC_DCTLAGtoHansen{1};
    xstar_DCTLAGtoHansen{i} = xstar{1};
    energyCost_nominal_DCTLAGtoHansen{i} = energyCost_nominal{1};
    ET_mgd_solutions_DCTLAGtoHansen{i} = ET_out{1};
    problemCaps{i} = probCaps{1};
    
end

save(filename)
    
end
