classdef elevation_profile
    % Description: Object that contains the physical information about a
    % single pipeline route.
    properties (GetAccess = public, SetAccess = private)
        units % indicates if we are using US or SI units (0 is SI and 1 is US)
        x % distance from starting point
        z % elevation at x points
        delta_x % change in x
        delta_z % change in z
        delta_L % change in length (Euclidean distance)
        ID % ID number of elevation profile (corresponds to GIS)
    end
    methods
        function obj = elevation_profile(L)
            % create an array of L pipelines
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj(1:L) = elevation_profile; % Preallocate object array
            end
        end
        function obj = set(obj, prop, value)
            % Set the properties of this instance of class
            if isequal(prop,'x')
                obj.x = value;
            elseif isequal(prop,'z')
                obj.z = value;
            elseif isequal(prop,'delta_x')
                obj.delta_x = value;    
            elseif isequal(prop,'delta_z')
                obj.delta_z = value;
            elseif isequal(prop,'delta_L')
                obj.delta_L = value;
            elseif isequal(prop,'ID')
                obj.ID = value;
            elseif isequal(prop,'units')
                obj.units = value;                
            else disp(['Property "' prop, '" not defined'])
            end
        end
    end
end