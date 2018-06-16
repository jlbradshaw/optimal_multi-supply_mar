% Copyright Jonathan L. Bradshaw 2018

classdef pumpstation_collection
    properties (GetAccess = public, SetAccess = private)
        construction_cost % construction cost ($)
        om_cost % Operation & maintenance cost ($/yr)
        energy_use % total system energy use (kWh/day)
        pumpstations % pumpstations that are part of the collection
        prvs % prvs that are part of the collection
        LCC % life cycle cost of the pumpstation_collection + pipes
        LCC_pumpstations % life cycle cost of the pumpstations
        LCC_pipes % life cycle cost of the pipes
    end
    methods
        function obj = pumpstation_collection(varargin)
            % create an array of L operation scenarios
            if nargin == 1 % Allow nargin == 0 syntax
                L = varargin{1};
                obj(1:L) = pumpstation_collection; % Preallocate object array
            elseif nargin == 2
                m = varargin{1};
                n = varargin{2};
                obj(1:m,1:n) = pumpstation_collection;
            end
        end
        function obj = set(obj, prop, value)
            % Set the properties of this instance of class
            % Can call either by obj = obj.set(prop, value) or obj = set(obj, prop, value)            
            if isequal(prop,'construction_cost')
                obj.construction_cost = value;
            elseif isequal(prop,'om_cost')
                obj.om_cost = value;
            elseif isequal(prop,'energy_use')
                obj.energy_use = value;
            elseif isequal(prop,'pumpstations')
                obj.pumpstations = value;
            elseif isequal(prop,'prvs')
                obj.prvs = value;
            elseif isequal(prop,'LCC')
                obj.LCC = value;
            elseif isequal(prop,'LCC_pumpstations')
                obj.LCC_pumpstations = value;
            elseif isequal(prop,'LCC_pipes')
                obj.LCC_pipes = value;
            else disp(['Property "' prop '" not defined in class pumpstation_collection'])
            end
        end
        function obj = compute_costs(obj) % Compute the total costs over an array of pumpstations (i.e., a pumpstation_collection)
            pstations = obj.pumpstations;
            c_cost = zeros(size(pstations));
            o_cost = zeros(size(pstations));
            e_use = zeros(size(pstations));
            for i = 1:length(pstations(1,:)) % loop through the cell array to find the cost associated with each pumpstation vector in the cell array of pumpstation_collection
                for j = 1:length(pstations(:,1))
                    temp = pstations{i,j};
                    c_cost(i,j) = sum([temp.construction_cost]);            
                    o_cost(i,j) = sum([temp.om_cost]);
                    e_use(i,j) = sum([temp.energy_use]);
                end
            end
            obj = set(obj,'construction_cost',c_cost);
            obj = set(obj,'om_cost',o_cost);
            obj = set(obj,'energy_use',e_use);
        end
        function obj = compute_pumping(obj,pipeline_x, pipe_op, psimax, psimin) % Compute the pumping and PRV requirements for a particular pipeline (i.e., 'pipeline_x' given the maximum and minumum pressure heads)                    
            ev = pipeline_x.elev_profile;            
            hfs = pipe_op.hfs_perL;
            pumpstation_cellArray = cell(size(hfs)); % initialize empty cell array to store the arrays of pumps and
            % prvs for each hfs option
            prv_cellArray = cell(size(hfs));
            % Iterate through each option of hfs, which depends on flowrate and pipe type (i.e.,
            % diameter and C)
            for i = 1:length(hfs(1,:))
                for j = 1:length(hfs(:,1))                    
                    [pump_needed, PRV_needed] = pumpdesign(ev, hfs(i,j), psimax, psimin);
                    
                    % Iterate through all the needed pumps and add them as
                    % pumpstation objects
%                     if isequal(pump_needed,{'NA'})
                    if ~isempty(pump_needed)% verify that pumps are needed
                    num_pumps = length(pump_needed(:,1));
                    pumps = pumpstation(num_pumps);
                    Q_afd = pipe_op.flowrates(j);
                    Q_gpm = units_check(units_check(Q_afd,'afd2cfs'),'cfs2gpm'); % maximum flow rate (gpm)
                    c_cost = 3.12*10^(0.7583*log10(Q_gpm)+3.1951);
                    om_cost = 10000 + 0.05*c_cost;
                    Q_afy_avg = units_check(Q_afd,'afd2afy'); % annual average flow in AFY
                    for k = 1:num_pumps
                        H_max = pump_needed(k,1);
                        H_ave = H_max; % use this assumption for now, might come back later to change
                        hp_ave = Q_afy_avg/1.613*H_ave/3956/0.75; % average brake horsepower
                        e_use = units_check(hp_ave*24,'hp2kw'); % daily energy use in kWh
                        pumps(k) = set(pumps(k),'head_max',pump_needed(k,1));
                        pumps(k) = set(pumps(k),'location_index',pump_needed(k,2));
                        pumps(k) = set(pumps(k),'construction_cost',c_cost);
                        pumps(k) = set(pumps(k),'om_cost',om_cost);
                        pumps(k) = set(pumps(k),'energy_use',e_use);
                        pumps(k) = set(pumps(k),'Q_max',Q_gpm);
                    end
                    pumpstation_cellArray{i,j} = pumps;
                    else pumpstation_cellArray{i,j} = 'NA';
                    end
                    
                    % Iterate through all the needed PRVs and add them as
                    % PRV objects
                    if ~isempty(PRV_needed)% verify that PRVs are needed
                        num_prvs = length(PRV_needed(:,1));
                        ps = prv(num_prvs);
                        for k = 1:num_prvs
                            ps(k) = set(ps(k),'head_max',PRV_needed(k,1));
                            ps(k) = set(ps(k),'location_index',PRV_needed(k,2));
                        end
                        prv_cellArray{i,j} = ps;
                    else prv_cellArray{i,j} = 'NA';
                    end
                end
                
            end
            % Assign the cell arrays to the pumpstation collection object                
            obj = set(obj,'pumpstations',pumpstation_cellArray);
            obj = set(obj,'prvs',prv_cellArray);
        end            
    end
end