% Copyright Jonathan L. Bradshaw 2018

classdef prv
    properties (GetAccess = public, SetAccess = private)
        construction_cost % construction cost ($)
        om_cost % Operation & maintenance cost ($/yr)
        energy_use % energy use (kWh/day)
        hp_ave % average horsepower of pump (horsepower)
        head_max %
        head_ave %
        location_index
    end
    methods
        function obj = prv(varargin)
            % create an array of L operation scenarios
            if nargin == 1 % Allow nargin == 0 syntax
                L = varargin{1};
                obj(1:L) = prv; % Preallocate object array
            elseif nargin == 2
                m = varargin{1};
                n = varargin{2};
                obj(1:m,1:n) = prv;
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
            elseif isequal(prop,'hp_ave')
                obj.hp_ave = value;
            elseif isequal(prop,'head_ave')
                obj.head_ave = value;
            elseif isequal(prop,'head_max')
                obj.head_max = value;
            elseif isequal(prop,'location_index')
                obj.location_index = value;
            end
        end
    end
end