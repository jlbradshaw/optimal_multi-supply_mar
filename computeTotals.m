% Copyright Jonathan L. Bradshaw 2018

% Compute ET, rw, and infil, etc. totals


color_bank15 = 1/255*[0,0,0; ... %black
    0,73,73; ... % dark green
    0,146,146; ... % dark blue-green
    255,109,182; ... % dark pink
    255,182,119; ... % light pink
    73,0,146; ... % dark purple
    0,109,219; ... % medium blue
    182,109,255; ... % lavender purple
    109,182,255; ... % sky blue
    182,219,255; ... % pastel blue
    146,0,0; ... % red-brown
    146,73,0; ... % brown
    219,109,0; ... % orange-brown
    36,255,36; ... % fluorescent green
    %     0,255,128; ... % softer fluorescent green
    % 0,255,0; ... % another green
    255,255,109]; % yellow

load('c:\Code and input\Results\final results\workspace_HTP_3-105mgd_30yr_fixedUncap-0_storage-1_fixedRWCap-1_ETinConst-0_SF-1_2018-01-28_12-45.mat')
ET_perYear = zeros(1,length(ET_mgd_solutions_HTPtoHansen));

for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    ET_perYear(i) = sum(ET_mgd_solutions_HTPtoHansen{i});
end

rw_total = zeros(1,length(ET_mgd_solutions_HTPtoHansen));
for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    rw_total(i) = sum(xstar_array_HTPtoHansen{i}(rw_indices));
end

infil_total = zeros(1,length(ET_mgd_solutions_HTPtoHansen));
for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    infil_total(i) = sum(xstar_array_HTPtoHansen{i}(infiltration_indices));
end

plot(ET_perYear./infil_total)

plot(rw_total./infil_total)


%% Compute number of underused WRF days
underusedWRF_days = zeros(1,length(ET_mgd_solutions_HTPtoHansen));
temp_threshold = 0.90;
for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    underusedWRF_days(i) = sum(xstar_array_HTPtoHansen{i}(rw_indices)<temp_threshold*max(xstar_array_HTPtoHansen{i}(rw_indices)));
end
plot(underusedWRF_days(6:end)/number_timesteps,'.')
temp_xlim = get(gca,'XLim')
set(gca,'XLim',[0, temp_xlim(2)])
temp_ylim = get(gca,'YLim')
set(gca,'YLim',[0, 1])
xlabel('WRF capacity (MGD) - 7')
ylabel(strcat('Fraction of days operating at <',num2str(temp_threshold*100),'% capacity'))

%% Compute fraction of full capacity used
fig2 = figure();
ax = axes('Parent',fig2);
% ax_compiled.Color = [1,1,1]*0.8;
set(ax,'nextplot','add');

set(gcf,'PaperPositionMode','auto'); % WYSIWYG
set(gca,'ygrid','on')
set(gca,'GridColor',[1,1,1]*0.8)
set(gca,'GridAlpha',0.5)

set(ax,'FontName','Helvetica')
set(ax,'FontSize',10)

marker_size = 4;
marker_edgewidth = 0.1;
load('c:\Code and input\Results\final results\workspace_HTP_3-105mgd_30yr_fixedUncap-0_storage-1_fixedRWCap-1_ETinConst-0_SF-1_2018-01-28_12-45.mat')
capacityUsed = zeros(1,length(ET_mgd_solutions_HTPtoHansen));
for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    capacityUsed(i) = sum(xstar_array_HTPtoHansen{i}(rw_indices))/number_timesteps/max(xstar_array_HTPtoHansen{i}(rw_indices));
end
temp_x = 8:97;
plot(ax,temp_x,capacityUsed(6:6-1+length(temp_x))*100,'ok','MarkerSize',marker_size,'LineWidth',marker_edgewidth,'Color',color_bank15(9,:))
xlabel('WRF capacity (MGD)')
ylabel('WRF utilization (%)')
set(gca,'YLim',[0,100])

load('c:\Code and input\Results\final results\workspace_HTP_3-105mgd_30yr_fixedUncap-1_storage-1_fixedRWCap-1_ETinConst-0_SF-1_2018-01-24_23-32.mat')

for i = 1:length(ET_mgd_solutions_HTPtoHansen)
    capacityUsed(i) = sum(xstar_array_HTPtoHansen{i}(rw_indices))/number_timesteps/max(xstar_array_HTPtoHansen{i}(rw_indices));
end
temp_x = 8:79;
plot(ax,temp_x,capacityUsed(6:6-1+length(temp_x))*100,'^k','MarkerSize',marker_size,'LineWidth',marker_edgewidth,'Color',color_bank15(12,:))
legend('Flexible potential to receive recycled water','Conservative potential to receive recycled water ','Location','southeast')


%% Energy
plot([energyCost_nominal_HTPtoHansen{1:end}],'.')