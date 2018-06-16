% Copyright Jonathan L. Bradshaw 2018

classdef pipe_op
    properties (GetAccess = public, SetAccess = private)
        length % the pipe_op
        units % indicates if we are using US (units = 1) or SI (units = 0) units
        materials % name of the material
        hfs_perL % frictional head loss in the pipe per unit length in ft or m. As defined in pipeline_routing, hfs is a 2-D array in which each column is a different flow rate and each row is a different pipe size and/or material.
        Cs % Hazen-Williams coefficient, which reflects the material
        diameters % inner diameter of the pipe in ft or m
        flowrates % flowrate in AFD
        energy_reqs % energy required to move the water as head (ft or m water)
        unit_cost_const % pipe cost to install and O&M per unit length
        unit_cost_OM % pipe cost O&M per unit length
        pumpstations % what pumpstations are needed
        prvs % what pressure reducing valves are needed
    end
    methods
        function obj = pipe_op(L)
            % create an array of L operation scenarios
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj(1:L) = pipe_op; % Preallocate object array
            end
        end
        function obj = set(obj, prop, value)
            % Set the properties of this instance of class
            % Can call either by obj = obj.set(prop, value) or obj = set(obj, prop, value)            
            if isequal(prop,'units')
                obj.units = value;
            elseif isequal(prop,'length')
                obj.length = value;
            elseif isequal(prop,'materials')
                obj.materials = value;
            elseif isequal(prop,'hfs_perL')
                obj.hfs_perL = value;
            elseif isequal(prop,'Cs')
                obj.Cs = value;
            elseif isequal(prop,'diameters')
                obj.diameters = value;
            elseif isequal(prop,'flowrates')
                obj.flowrates = value;
            elseif isequal(prop,'energy_reqs')
                obj.energy_reqs = value;
            elseif isequal(prop,'unit_cost_const')
                obj.unit_cost_const = value;
            elseif isequal(prop,'unit_cost_OM')
                obj.unit_cost_OM = value;
            elseif isequal(prop,'pumpstations')
                obj.pumpstations = value;
            elseif isequal(prop,'prvs')
                obj.prvs = value;
            else disp(['Property "' prop '" not defined in class pipe_op'])
            end
        end        
        
        function obj = import_pipe_types(obj, filename)
            fileID = fopen(filename);
            C = textscan(fileID, '%s %u %f %f %f %f','Delimiter',',','HeaderLines',1); % Assumes that there is one header line
            fclose(fileID);
         
            obj = set(obj,'materials', C{1}); % assumes first row is headers, so the text data are offset by 1 row. textdata row must be offset by 1. The data rows (i.e., numerical data) do not need to be offset.
            obj = set(obj,'units', C{2}); % units is 0 for SI and 1 for US
            
            obj = set(obj,'Cs', C{3});
            obj = set(obj,'diameters', C{4});
            obj = set(obj, 'unit_cost_const', C{5});
            obj = set(obj, 'unit_cost_OM', C{6});

        end           
    end
end