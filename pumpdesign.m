% Copyright Jonathan L. Bradshaw 2018

function [pipeLength_ft, pump_needed_no0, PRV_needed_no0] = pumpdesign(elev_profile_m, hfs_perL_ft, psiMax_ft, psiMin_ft)

% Purpose: Design a pump system that minimizes the number and size of pumps (i.e., minimizes costs) and meets the
% required performance criteria (i.e., pressure between psimin and psimax (e.g., 200 psi), and
% there is sufficient head convey the water from recycled water facility to
% spreading basin.

% Inputs: Elevation profile (in m), frictional headloss hf (foot head loss per foot pipe length), maximum allowable pressure head psimax (in ft),  minimum allowable pressure head psi min (in ft)

% Outputs:
% (0) Length of the pipeline in feet
% (1) Array that contains the pump head that is needed and the location along pipeline
% (2) Array that contains the Pressure Release Valves head that are needed and
% the location along the pipeline

len_convert = units_check(0,1,'length'); % convert elevation profile from m to ft
pump_needed = []; % empty matrix where first column indicates size of pump needed and second column is distance along pipeline that pump is needed
PRV_needed = []; % similar to 'pumps_needed' matrix; empty matrix where first column indicates head dissipation requirement of pressure-reducing valve (PRV) and second column is distance along pipeline that the PRV is needed


x_ft = elev_profile_m(:,1)*len_convert; % Cumulative distance from start in feet.
z_ft = elev_profile_m(:,2)*len_convert; % Elevation in feet

delta_x_ft = x_ft(2:end) - x_ft(1:end-1);
delta_z_ft = z_ft(2:end) - z_ft(1:end-1);
delta_L_ft = delta_x_ft;
pipeLength_ft = sum(delta_L_ft);

psi_ft = zeros(size(z_ft)); % pressure head vector - defined starting at the spreading basin and ending at the recycled water facility
psi_ft(1) = 0;

% Working forwards from the WRF to SB
% Bernoulli Equation for constant velocity: z1 + psi1 = z2 + psi2 + h_L
for i = 1:(length(z_ft)-1)
    z1_ft = z_ft(i);
    z2_ft = z_ft(i+1);
    psi1_ft = psi_ft(i);
    h_L_ft = hfs_perL_ft*(delta_L_ft(i)); % Compute head loss due to friction in the segment
    psi2 = (z1_ft - z2_ft) + psi1_ft - h_L_ft; % rearrangement of Bernouilli equation to solve for psi2
    
    %     confirm that psi2 is within permissible range
    if psi2 >= psiMax_ft % if psi2 exceeds allowable pressure      
        PRV_needed = [PRV_needed; psi2-psiMax_ft, i]; % add a PRV - might figure out a better way to perform this rather than by dynamic allocation
        psi1_ft = psi1_ft - (psi2-psiMax_ft);
    elseif psi2 <= psiMin_ft % if psi2 pressure is too low, we must install another pump
        if i == 1
            pump_needed = [pump_needed; psiMax_ft-(psi1_ft), i];
        else
            pump_needed = [pump_needed; psiMax_ft-(psi1_ft), i]; % add a pump - might figure out a better way to perform this rather than by dynamic allocation           
        end
        psi1_ft = psiMax_ft; % install a pump to increase pressure head to maximum permissible pressure
    end
    psi_ft(i) = psi1_ft;
    psi2 = (z1_ft - z2_ft) + psi1_ft - h_L_ft; % rearrangement of Bernouilli equation to solve for psi2
    psi_ft(i+1) = psi2;
    
end

if ~isempty(pump_needed)
    
    if ~isempty(PRV_needed)
        pump_segs_index = pump_needed(:,2); % index location of pumps
        PRV_segs_index = PRV_needed(:,2); % index location of PRVs
        
        for j = 1:length(PRV_segs_index) % for each location of the PRV, check the upstream segment to see how much pressure could be reduced, and subtract as much head as possible (e.g., subtract the PRV head reduction from the last pump before it IFF there is excess pressure head (i.e., make sure that psi > psimin for all segments between the PRV and pump)
            PRV = PRV_needed(j,:); % the PRV we are looking at right now
            pump_index = find(pump_segs_index<PRV_segs_index(j),1,'last'); % find the last pump upstream of the PRV
            
            % isolate the pressure head section between the PRV and the upstream pump
            section = pump_needed(pump_index,2):PRV(2);
            psi_section = psi_ft(section);
            
            % determine the excess pressure head = the minimum of the
            % a.) the difference between the minimum pressure head in
            % the segment and psimin and b.) the head dissipated by the
            % PRV [PRV(1)]
            excess_psi = min((min(psi_section)- psiMin_ft),PRV(1));
            
            % change the pump head and PRV head required; also change the
            % pressure heads in the affected section
            if (pump_needed(pump_index,1) - excess_psi) > 0
                pump_needed(pump_index,1) = pump_needed(pump_index,1) - excess_psi;
                PRV_needed(j,1) = PRV(1) - excess_psi;
                psi_ft(section) = psi_section - excess_psi;
            end
                       
        end
    end
    
    last_pump_index = pump_needed(end,2); % find the index location of the last pump needed
    
    nomorePRVs = 1;
    if ~isempty(PRV_needed) % if there are PRVs
        nomorePRVs = 0;
        if last_pump_index > PRV_needed(end,2) % if the last pump is after the last PRV
            nomorePRVs = 1;
        elseif sum(PRV_needed((last_pump_index:end),1)) == 0 % or if there were PRVs, but they are all non-operational
            nomorePRVs = 1;
        end
    end
        
    if nomorePRVs
        % If there are no PRVs after the last pump, resize the last pump so we don't have excess pressure head at the
        % outlet
        % isolate the pressure head section between the PRV and the upstream pump
        psi_section = psi_ft(pump_needed(end,2):end);
        min_psi = min(psi_section);
        if (sum(min_psi == psi_section) == 1) && (find(psi_section == min_psi) == length(psi_section))% if there is only one place where the pressure is a minimun and it is the outlet (i.e., the last element in the psi_section vector), then set outlet pressure head to zero
            last_pump_head = pump_needed(end,1) - min_psi;
        else last_pump_head = pump_needed(end,1) - (min_psi - psiMin_ft); % calculate the head needed for the last pump
        end
        pump_needed(end,1) = last_pump_head; % change the value of the last pump needed to calculated value
        
        % change the pressure head computations given the new pressure head
        % values
        psi_ft(last_pump_index) = last_pump_head; % change the value of the pressure head at the location of the last pump head
        for j = last_pump_index:(length(z_ft)-1)
            z1_ft = z_ft(j);
            z2_ft = z_ft(j+1);
            psi1_ft = psi_ft(j);
            h_L_ft = hfs_perL_ft*(delta_L_ft(j)); % convert head loss due to friction (assumed ft) to m
            psi2 = (z1_ft - z2_ft) + psi1_ft - h_L_ft; % rearrangement of Bernouilli equation to solve for psi2
            psi_ft(j+1) = psi2;
        end
    end
    
    % remove any pumps or PRVs with 0 head.
    if isempty(pump_needed)
        pump_needed_no0 = pump_needed;
    else pump_needed_no0 = pump_needed(pump_needed(:,1)~=0,:);
    end
    if isempty(PRV_needed)
        PRV_needed_no0 = PRV_needed;
    else
        PRV_needed_no0 = PRV_needed(PRV_needed(:,1)~=0,:);
    end       
end