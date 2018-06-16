% Copyright Jonathan L. Bradshaw 2018

function output = units_check(varargin)
% Purpose: Check to confirm that units are consistent (nargin ==3). Convert
% units (nargin == 2). Return the number of days per month in water year
% order (nargin == 1).
if nargin == 3
%     format: varargin = {{b},{a},{param}}. 'output' is a conversion factor
%     to multiply a by to get units of a. (e.g., output converts [a] to
%     [b])
        
    b = varargin{1};
    a = varargin{2};    
    param = varargin{3};
    if isequal(param,'length')
        c = 3.28084; % unit conversion for meters to feet (i.e., 1 m = c ft)
        if a == b % if a and b are the same, then conversion is not necessary (i.e., conversion_factor is 1
            conversion_factor = 1;
        elseif a == 0  % else, if a is SI and b is US
            conversion_factor = 1/c;
        elseif a == 1
            conversion_factor = c;
        else disp('Units error')
        end
    end
    output = conversion_factor;
elseif nargin == 2
    a = varargin{1};
    b = varargin{2};
    if isequal(b,'afm2cfs') % convert a from acre-feet per month to cubic-feet per second
        c = 60.3706114;
        output = a/c;
    elseif isequal(b,'cfs2afm')
        c = 60.3706114;
        output = a*c;
    elseif isequal(b,'afd2cfs')
        c = 1.98347107;
        output = a/c;
    elseif isequal(b,'afd2afy')
        c = 365.25;
        output = a*c;
    elseif isequal(b,'afy2afd')
        c = 365.25;
        output = a/c;
    elseif isequal(b,'cfs2afd')
        c = 1.98347107;
        output = a*c;
    elseif isequal(b,'afy2cfs')
        c = 12*60.3706114;
        output = a/c;
    elseif isequal(b,'cfs2afy')
        c = 12*60.3706114;
        output = a*c;
    elseif isequal(b,'cfs2gpm')
        c = 448.831169;
        output = a*c;
    elseif isequal(b,'gpm2cfs')
        c = 448.831169;
        output = a/c;
    elseif isequal(b,'afy2mgd')
        c = 1120.88568;
        output = a/c;
    elseif isequal(b,'mgd2afy')
        c = 1120.88568;
        output = a*c;
    elseif isequal(b,'afm2afy')
        c = 12;
        output = a*c;
    elseif isequal(b,'afy2afm')
        c = 12;
        output = a/c;
    elseif isequal(b,'afd2mgd')
        c = 3.06888328;
        output = a/c;
    elseif isequal(b,'mgd2afd')
        c = 3.06888328;
        output = a*c;
    elseif isequal(b,'hp2kw')
        c = 1.34102209;
        output = a/c;
    elseif isequal(b,'kw2hp')
        c = 1.34102209;
        output = a*c;
    elseif isequal(b,'mg2g')
        c = 1e6;
        output = a*c;
    elseif isequal(b,'g2mg')
        c = 1e6;
        output = a/c;
    elseif isequal(b,'gpm2afy')
        c = 1.61301;
        output = a*c;
    elseif isequal(b,'mgd2gpm')
        c = 1e6/(60*24);
        output = a*c;
    elseif isequal(b,'afy2gpm')
        c = 1.61301;
        output = a/c;
    elseif isequal(b,'afy2mly')
        c = 1.23348;
        output = a*c;
    elseif isequal(b,'kwh2j')
        c = 3.6e6;
        output = a*c;
    elseif isequal(b,'j2kwh')
        c = 3.6e6;
        output = a/c;
    elseif isequal(b,'af2m3')
        c = 1233.48;
        output = a*c;
    elseif isequal(b,'month2days')
        c = 30;
        output = a*c;
    elseif isequal(b,'year2month')
        c = 12;
        output = a*c;
    elseif isequal(b,'month2year')
        c = 12;
        output = a/c;
    elseif isequal(b,'year2days')
        c =365.25;
        output = a*c;
    elseif isequal(b,'s2d')
        c = 1.15741e-5;
        output = a*c;
    elseif isequal(b,'d2s')
        c = 1.15741e-5;
        output = a/c;
    elseif isequal(b,'m3s2mgd')
        c = 22.82;
        output = a*c;
    elseif isequal(b,'mgd2m3s')
        c = 22.82;
        output = a/c;
    elseif isequal(b,'ft2m')
        c = 0.3048;
        output = a*c;
    elseif isequal(b,'m2ft')
        c = 0.3048;
        output = a/c;
    elseif isequal(b,'af2mg')
        c = 0.325851;
        output = a*c;
    elseif isequal(b,'mg2af')
        c = 0.325851;
        output = a/c;
    elseif isequal(b,'afm2mgd')
        c = 10712.9/1e6;
        output = a*c;
    elseif isequal(b,'mgd2afm')
        c = 10712.9/1e6;
        output = a/c;
    elseif isequal(b,'cfs2mgd')
        c = 646317/1e6;
        output = a*c;
    elseif isequal(b,'mgd2cfs')
        c = 646317/1e6;
        output = a/c;
    elseif isequal(b,'mg2m3')
        c = 3785;
        output = a*c;
    elseif isequal(b,'m32mg')
        c = 3785;
        output = a/c;
    end
elseif nargin == 1
    output = [31,30,31,31,28.25,31,30,31,30,31,31,30];    
    
else disp('Incorrect number of inputs for "converstion_factor"') 
end