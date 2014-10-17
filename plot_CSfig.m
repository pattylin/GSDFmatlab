function plot_CSfig(event,CS,CSplot)

% Function to plot representative figure for CSmeasure script. CSplot is
% an arry that contains most of the information used to plot the figure
% Get information concnering station locations
load stalst;
issaved = 0;
% Pretty colors!
ORANGE=[255 69 0]./255;
BLUE=[44/255 77/255 143/255];
RED=[231/255 47/255 39/255];
GREEN=[19/255 166/255 50/255];

if ~exist('periods','var')
    setup_parameters;
    periods = parameters.periods;
    lalim = parameters.lalim;
    lolim = parameters.lolim;
end

sta1 = CS.sta1;
sta2 = CS.sta2;
sta1_stnm = event.stadata(sta1).stnm;
dist1 = event.stadata(sta1).dist;
sta2_stnm = event.stadata(sta2).stnm;
dist2 = event.stadata(sta2).dist;

% Get information from CSplot array
data1 = CSplot.data1;
win_data2 = CSplot.win_data;
taxis1 = CSplot.taxis1;
taxis2 = CSplot.taxis2;
nband_win_xcors = CSplot.nband_win_xcors;
lag = CSplot.lag;
xcor = CSplot.xcor;
win_xcor = CSplot.win_xcor;

%% Plot a Map of Station Locations
rayx = [event.stadata(sta1).stla;event.stadata(sta2).stla];
rayy = [event.stadata(sta1).stlo;event.stadata(sta2).stlo];
f46 = figure(46);
clf
subplot(4,2,[1,3])
ax = worldmap(lalim,lolim);
plotm(stla,stlo,'ok','markerfacecolor',GREEN);
plotm(event.stadata(sta1).stla,event.stadata(sta1).stlo,'ok','markerfacecolor',RED);
plotm(event.stadata(sta2).stla,event.stadata(sta2).stlo,'ok','markerfacecolor',RED);
plotm(rayx,rayy,'-k');
title('Station Pair')

%% Plot the Raw and Windowed Data
subplot(4,2,2)
plot(taxis1,data1);
xlim([0 dist2/2])
text(500,max(data1)-0.5*max(data1),sta1_stnm,'color',ORANGE,'fontweight','bold');
title('Station 1');

subplot(4,2,4)
plot(taxis2,win_data2);
xlim([0 dist2/2])
text(500,max(win_data2)-0.5*max(win_data2),sta2_stnm,'color',ORANGE,'fontweight','bold');
title('Station 2');

%% Plot the Power of the Cross Correlation Across Periods
subplot(4,2,[5,7])
[xi yi] = ndgrid(lag,periods);
for ip = 1:length(periods)
    norm_nbands(:,ip) = nband_win_xcors(:,ip)./max(abs(nband_win_xcors(:,ip)));
end
contourf(xi,yi,norm_nbands);
hold on
xlim([-3*max(periods) 3*max(periods)]);

goodind = find(CS.isgood == 1);
for ig = 1:length(goodind)
    plot(300,periods(goodind(ig)),'rp','markersize',15,'markerfacecolor',RED)
    %text(-250,periods(goodind(ig)),'Good','color',ORANGE,'fontweight','bold');
end
hold off
%% Plot the Cross Correlograms
subplot(4,2,6)
plot(lag,xcor);
xlim([-500 500])
title('Cros Corr.');

subplot(4,2,8)
plot(lag,win_xcor);
xlim([-500 500])
title('Windowed Cross Corr.');

pos = get(f46,'Position');
set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
    'Position',[pos(1) pos(2) 700 660]);
drawnow;

if issaved ==1
    FIG46_id = [event.id,'_',sta1_stnm,'_',sta2_stnm,'.jpg'];
    print(f46,'-djpeg',FIG46_id);
end
end
