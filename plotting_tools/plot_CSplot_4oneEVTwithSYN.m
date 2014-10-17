function plot_CSplot_4oneEVT(EVT)


% Function to plot representative figure for CSmeasure script. CSplot is
% an arry that contains most of the information used to plot the figure
% Get information concnering station locations
%mk_stalist
load stalst;
issaved = 1;
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
    comp = parameters.component;
	gridsize = parameters.gridsize;

	xnode=lalim(1):gridsize:lalim(2);
	ynode=lolim(1):gridsize:lolim(2);
	Nx=length(xnode);
	Ny=length(ynode);
	[xi yi]=ndgrid(xnode,ynode);
	
end

load(['./CSmeasure/',num2str(EVT),'_cs_',comp,'.mat'])
load(['./eventmat/',num2str(EVT),'_',comp,'.mat'])


CSnum = length(eventcs.CS); % CSnum = all pairs number for this events

%% Plot a Map of Station Locations for all pairs.
figure(45)
clf
N=3; M = floor(length(periods)/N)+1;
clf
for ip = 1:length(periods)
   subplot(M,N,ip)
   ax = worldmap(lalim, lolim);
   set(ax, 'Visible', 'off')
   if exist('drawlocal.m','file')
         drawlocal
   end
 
    temparray = vertcat(eventcs.CS.isgood);
 	goodind = find(temparray(:,ip) == 1 );
     %%
 
 	for ig = 1:length(goodind)
 		CS = eventcs.CS(ig);
 		CSplot = eventcs.CSplot(ig);
 		sta1 = CS.sta1;
 		sta2 = CS.sta2;
 		sta1_stnm = event.stadata(sta1).stnm;
 		sta2_stnm = event.stadata(sta2).stnm;
 	
 		rayx = [event.stadata(sta1).stla;event.stadata(sta2).stla];
 		rayy = [event.stadata(sta1).stlo;event.stadata(sta2).stlo];
 		plotm(stla,stlo,'ok','markerfacecolor',GREEN);
 		plotm(rayx,rayy,'-k');
		plotm(xi,yi,'vk')
     end
     title([num2str(length(goodind)) , 'Npairs in period  ',num2str(periods(ip))])
 
 end
 drawnow;





%for iCSnum = 4:4
for iCSnum = 1:CSnum
    
    iCSnum
    CS = eventcs.CS(iCSnum);
    CSplot = eventcs.CSplot(iCSnum);
	avgphv = eventcs.avgphv;

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
    [xi yi] = ndgrid(CSplot.lag,periods);

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
    plot(CSplot.lag,xcor);
    xlim([-500 500])
    title('Cros Corr.');

    subplot(4,2,8)
    plot(CSplot.lag,win_xcor);
    xlim([-500 500])
    title('Windowed Cross Corr.');

    pos = get(f46,'Position');
    set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 660]);
    drawnow;
    

    f47 = figure(47);
    clf
    for ip = 1:length(periods)
        subplot(length(periods),1,ip)
        if ( CS.isgood(ip) == 1 ) 
        	plot(CSplot.lag,nband_win_xcors(:,ip),'b-');
        else
        	plot(CSplot.lag,nband_win_xcors(:,ip),'-','color', [0.5 0.5 0.5]);
		end
        title(['CSisgood =  ', num2str(CS.isgood(ip)),'Period= ',num2str(periods(ip)),' (s)',' avgphv: ',num2str(avgphv(ip))]);
    end
    pos = get(f47,'Position');
    set(gcf, 'PaperOrientation', 'portrait','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 900]);
    drawnow;
        

    if issaved ==1
        FIG46_id = [event.id,'_',sta1_stnm,'_',sta2_stnm,'.jpg'];
        print(f46,'-djpeg',FIG46_id);
        FIG47_id = [event.id,'_',sta1_stnm,'_',sta2_stnm,'.ps'];
        print(f47,'-dpsc2',FIG47_id);
        
    end
end
end
