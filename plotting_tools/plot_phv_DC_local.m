% run after stack_helmholtz or stack_phv 
% Plot dispersion curve for average phase velocity
% Plot choice: 
%		click a point on the map, only plot the DC for that point
%               for each event and avg result
% pylin.patty 2013.11.24
%
% add for after ekional stacking 
% pylin.patty 2014.0421

clear all



SW_type = 'R'  % datatype: 'R': rayleigh;'L': Love

addpath('/Users/pappi/Research/PGS/MATLABdatalib')


% debug setting
isfiguretoPS = 1;
ishelmholtz = 0; % 1 more helmo, 0 for eikonal

if isfiguretoPS 
dispersion_curve_output_path = './DispersionCurve_local/';
if ~exist(dispersion_curve_output_path)
   mkdir(dispersion_curve_output_path)
end
end 


setup_parameters;
comp = parameters.component;
periods = parameters.periods;

lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;


xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);



% ===================================
if ishelmholtz
    %helmholtz_file = ['helmholtz_stack_',comp,'.mat']
    load(['helmholtz_stack_',comp,'.mat'])
else
    %eikonal_file = ['helmholtz_stack_',comp,'.mat']
    load(['eikonal_stack_',comp,'.mat'])
end


figure(11)
clf
title('stack for structure phv')
r = 0.03;
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(periods)
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    h1=surfacem(xi,yi,avgphv(ip).GV);
    %h1=surfacem(avgphv(ip).xi,avgphv(ip).yi,avgphv(ip).GV);
    % set(h1,'facecolor','interp');
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    avgv = nanmean(avgphv(ip).GV(:)); 
    %avgv = nanmean(avgphv(ip).GV(25))
    if isnan(avgv)
        continue;
    end
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    load seiscmap
    colormap(seiscmap)
end
drawnow;
[mlat mlon] = inputm(1);
[temp ilat] = min(abs(mlat - xnode));
[temp ilon] = min(abs(mlon - ynode));
plotm(mlat, mlon,'kx','linewidth',3,'markersize',10)
% ==================================================
% in stack_Helm save GV_cor_mat4plot.mat  GV_cor_mat
% in stack_phv save GV_matplot.mat GV_mat
if ishelmholtz
    load GV_cor_mat4plot.mat
    numbers_events = size(GV_cor_mat(:,:,:,:),3);
    GVnow = GV_cor_mat;
else
    load GV_mat4plot.mat
    numbers_events = size(GV_mat(:,:,:,:),3);
    GVnow = GV_mat;
end


load useEVT

for ie = 1:numbers_events
    figure(62)
    clf
    for ip=1:length(periods)
        evtavgphv(ip).phv(ie)= GVnow(ilat,ilon,ie,ip);
        plot(periods(ip),evtavgphv(ip).phv(ie),'.b','markersize',15);
        hold on;
    end    
    
    load pa5pvR
    plot(2*pi./w, phv,'r')

    xlim([10 110]);
    ylim([4.0 4.2]);
    xlabel('Period (s)');
    ylabel('Average Phase Velocity');
    gcarc(ie) = distance(eventinfo(ie).evla, eventinfo(ie).evlo, 9, -146);
    title(sprintf('Event %s dist %s',num2str(eventinfo(ie).id), num2str(gcarc(ie))));
    drawnow;
    if isfiguretoPS
       eventnamePS = [dispersion_curve_output_path,'/',num2str(eventinfo(ie).id),'.dispersion_curve.ps']
       %print( 'depsc', eventname); doesn't work 
       %print -dpsc2 eventnamePS    doesn't work
        print('-dpsc2',eventnamePS)
    end


end

for ip = 1:length(periods)
    sumphv(ip).phv = nanmean(evtavgphv(ip).phv);
    sumphv(ip).phvstd = nanstd(evtavgphv(ip).phv);
end


figure(63)
clf


plot(periods,[sumphv.phv],'xb','linewidth',2);
hold on
errorbar(periods,[sumphv.phv],[sumphv.phvstd],'xb', 'linewidth',1);

%load allgrid
%plot(periods,phvallgrid,'xg')

load pa5pvR
plot(2*pi./w, phv,'r')
hold on

load age_52_110Myr_iso_pvR
plot(period,phv,'k')

hold on;
load age_52_110Myr_aniso_pvR
plot(period,phv,'k:')



xlim([10 110]);
ylim([4.0 4.2]);
ylabel('Average Velocity (km/s)');
xlabel('Period (s)')
 if isfiguretoPS
    AVG_dispersion_curvePS = [dispersion_curve_output_path,'/AVG.dispersion_curve.ps']
    %print -dpsc2 'AVG_dispersion_curve.ps' work
     print('-dpsc2',AVG_dispersion_curvePS)
 end
ilat
ilon
