% Copyright Jonathan L. Bradshaw 2018

function elev_prof = import_piperoutes(filename)
% Input: Filename of piperoutes 'filename'
% Output: Cell Array of elevation profiles for pipe routes
% Dependencies: Internal MATLAB functions for importing data from a *.csv
% file
% Purpose: Read in a pipeline data from a *.csv file and assign it to a cell array

% read in data from CSV file
rawData = importdata(filename, ',');
x_index = 1;
z_index = 2;
lineID_index = 3;
line_ID = rawData.data(:,lineID_index); % column with route line ID numbers

% Need to separate each pipeline route
routes = unique(line_ID); % array of the unique routes numbers
elev_prof = cell(length(routes),2);

for i = 1:length(routes) % loop over all the routes
    elev_prof{i,1} = routes(i);
    rows = line_ID == routes(i); % isolate the indices of the rows that contain the desired route number
    
    x = rawData.data(rows,x_index); % cumulative distance
        
    z = rawData.data(rows,z_index); % change in elevation    
    
    elev_prof{i,2} = [x,z];
end