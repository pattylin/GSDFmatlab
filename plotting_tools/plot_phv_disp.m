% Plot dispersion curve for average phase velocity
% original got from Natalie 2013
% quick esitamte dispersion curve fron CSmeasure not from Eikonal
% pylin.patty 11/24,2013

 
clear all

% debug setting
isfiguretoPS = 1;

if isfiguretoPS 
dispersion_curve_output_path = './DispersionCurve_fromCSmeasure/';
if ~exist(dispersion_curve_output_path)
   mkdir(dispersion_curve_output_path)
end
end 



setup_parameters;
comp = parameters.component;
periods = parameters.periods;

%% Get information about only the used events
uevents = 'ueventlist.txt';
fid = fopen(uevents,'r');
if (fid == -1)
    error (['    Cannot open file: ', uevents]);
end
f = textscan(fid,'%f %s');
fclose(fid);

temp = f{2};
count = 0;
for iu=1:(length(temp)-1)
    if isequal(temp(iu),temp(iu+1)) == 1
        continue
    else
        count=count+1;
        %         disp(temp(iu));
        uevent(count) = temp(iu);
    end
end

%% Start to plot dispersion
% input path
eventcs_path = './CSmeasure/';
numbers_events = length(uevent)
for ie = 1:numbers_events
    csmatfiles(ie) = dir([eventcs_path,'/',char(uevent(ie)),'_cs_',comp,'.mat']);
    clear eventphv
    % read in data and set up useful variables
    temp = load([eventcs_path,csmatfiles(ie).name]);
    eventcs =  temp.eventcs;
    %     if exist('eventcs.avggrv','var') == 0
    %         continue;
    %     end
    disp(eventcs.id)
    eventinfo(ie).evla = eventcs.evla;
    eventinfo(ie).evlo = eventcs.evlo;
    eventinfo(ie).evdp = eventcs.evdp;
    
    figure(62)
    clf
    for ip = 1:length(periods)
        evtavgphv(ip).phv(ie) = eventcs.avgphv(ip); 
        plot(periods(ip),evtavgphv(ip).phv(ie),'.b','markersize',15);
        hold on
    end

    
    load pa5phv
    plot(2*pi./w, phv,'r')

    xlim([10 110]);
    ylim([4.0 4.2]); 
    xlabel('Period (s)');
    ylabel('Average Phase Velocity');
    title(sprintf('Event %s dist %s',num2str(eventcs.id), num2str(km2deg(mean(eventcs.dists))) ));
    drawnow;

    if isfiguretoPS
       eventnamePS = [dispersion_curve_output_path,'/',num2str(eventcs.id),'.dispersion_curve.ps']
       %print( 'depsc', eventname); doesn't work 
       %print -dpsc2 eventnamePS    doesn't work
        print('-dpsc2',eventnamePS)
    end
end

for ip = 1:length(periods)
    sumphv(ip).phv = nanmean(evtavgphv(ip).phv);
    sumphv(ip).phvstd = nanstd(evtavgphv(ip).phv);
end


% Look at the total average dispersion curves, but maybe this is not
% the right thing to look at as each event will vary slightly depending on
% the stations used
% 
 figure(63)
 clf

 plot(periods,[sumphv.phv],'xb','linewidth',2);
 hold on
 errorbar(periods,[sumphv.phv],[sumphv.phvstd],'xb', 'linewidth',1);



 load pa5phv
 plot(2*pi./w, phv,'r')
 xlim([10 110]);
 ylim([4.0 4.2]);
 ylabel('Average Velocity (km/s)');
 xlabel('Period (s)');

 if isfiguretoPS
    AVG_dispersion_curvePS = [dispersion_curve_output_path,'/AVG.dispersion_curve.ps']
    %print -dpsc2 'AVG_dispersion_curve.ps' work
     print('-dpsc2',AVG_dispersion_curvePS)
 end
