% Copyright Jonathan L. Bradshaw 2018

function [] = runGeneralScenarios_paper2

clear all

%% Select a scenario
scenario = 'DCT';
% scenario = 'LAG';
% scenario = 'HTP';
% scenario = 'LAG-DCT_series';
% scenario = 'LAG-DCT_parallel';
assessmentPeriod_yr = 30;
% assessmentPeriod_yr = 1;

is_uncap_fixed = false; % True if we fix spreading basin unused capacity (e.g., conservative paradigm); False if we want flexible paradigm
is_storage = true;                                                              % True if we want to consider ponding storage capacity
is_minRWCapacity_fixed = true; % Choose to fix either the RW production capacity (true) or the minimum amount of produced water (false)
is_ETinConstraints = false; % Choose to put ET in constraints
spreadingBasin_safetyFactor = 1; % Safety Factor that applies to spreading basin storage and infiltration rates

number_parpoolWorkers = 6;

    if strcmp('DCT', scenario)
        startCapacity_mgd = 3;                 
        rw_capacitymax = 27;
        runSingleWRFScenario;
    elseif strcmp('LAG', scenario)
        startCapacity_mgd = 3;
        rw_capacitymax = 40;
        runSingleWRFScenario;
    elseif strcmp('HTP', scenario)
        startCapacity_mgd = 3;        
        rw_capacitymax = 105;
        runSingleWRFScenario;
    elseif strcmp('LAG-DCT_series', scenario)
        startCapacity_mgd = [27,3];
        rw_capacitymax = [27,40];
        rw_capacitymax = [27,35];
        runTwoWRFScenario;
    elseif strcmp('LAG-DCT_parallel', scenario)
        startCapacity_mgd = [27,3];
        rw_capacitymax = [27,40];
        runTwoWRFScenario;
    end    


    function [] = runSingleWRFScenario()
                
        rw_capacityVector_mgd = startCapacity_mgd:rw_capacitymax;        
        
        
        rw_capacityVector_mgd = rw_capacityVector_mgd(rw_capacityVector_mgd <= rw_capacitymax); % need this for the last run to make sure we don't exceed the maximum capacity
        
        filename = strcat('results\workspace_',scenario,'_',num2str(rw_capacityVector_mgd(1)),'-',num2str(rw_capacityVector_mgd(end)),'mgd_',num2str(assessmentPeriod_yr),'yr_fixedUncap-',num2str(is_uncap_fixed),'_storage-',num2str(is_storage),'_fixedRWCap-',num2str(is_minRWCapacity_fixed),'_ETinConst-',num2str(is_ETinConstraints),'_SF-',num2str(spreadingBasin_safetyFactor),'_',datestr(now,'yyyy-mm-dd_HH-MM'))
        
        generalScenarios_paper2_par(rw_capacityVector_mgd, assessmentPeriod_yr, scenario, is_uncap_fixed, is_minRWCapacity_fixed, is_storage, is_ETinConstraints, spreadingBasin_safetyFactor, filename)
        
    end

    function [] = runTwoWRFScenario
        
        % The looping through capacities will be perfomed in the opt_paper
        % function, so pass through the entire rw_capacityVectors.
        rw_capacityVector_mgd = cell(1,2);
        rw_capacityVector_mgd{1} = startCapacity_mgd(1):rw_capacitymax(1);
        rw_capacityVector_mgd{2} = startCapacity_mgd(2):rw_capacitymax(2);
                
        filename = strcat('results\workspace_',scenario,'_upstream',num2str(rw_capacityVector_mgd{2}(1)),'-',num2str(rw_capacityVector_mgd{2}(end)),'mgd_terminal',num2str(rw_capacityVector_mgd{1}(1)),'-',num2str(rw_capacityVector_mgd{1}(end)),'mgd_',num2str(assessmentPeriod_yr),'yr_fixedUncap-',num2str(is_uncap_fixed),'_storage-',num2str(is_storage),'_fixedRWCap-',num2str(is_minRWCapacity_fixed),'_ETinConst-',num2str(is_ETinConstraints),'_SF-',num2str(spreadingBasin_safetyFactor),'_',datestr(now,'yyyy-mm-dd_HH-MM'))
                
        generalScenarios_paper2_par(rw_capacityVector_mgd, assessmentPeriod_yr, scenario, is_uncap_fixed, is_minRWCapacity_fixed, is_storage, is_ETinConstraints, spreadingBasin_safetyFactor, filename)
               
        
    end

end
