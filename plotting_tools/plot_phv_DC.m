%plot_phv_DC
% Plot dispersion curve for average phase velocity
% Plot choice: 
% 	isRegion_weighted: input the range of lat, lon,
%			   only calculate avg from these grids	 
%       isSTD_weighted   : only calcylate the avg when std <= 0.2
%       can work on both types
% pylin.patty 2013.11.24

clear all


SW_type = 'R'
% datatype:
% 'R': rayleigh
% 'L': Love

addpath('/Users/pappi/Research/PGS/MATLABdatalib')

% debug setting
isfiguretoPS = 1;
isRegion_weighted = 1; 
isSTD_weighted = 0;
ishelmholtz = 0; % 1 more helmo, 0 for eikonal 

if isfiguretoPS 
    dispersion_curve_output_path = './DispersionCurve_grid4/';
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





for ip=1:length(periods)
    GV(ip).wighted_region = 1;
    if isRegion_weighted 
        %grid 25 1by1;grid=4 2by 2
        GV(ip).wighted_region(1:3,1:3) = 0
		GV(ip).wighted_region(2,2)=1  
	    %gird 25	
        %GV(ip).wighted_region = xi < 10  & xi >= 9 & yi <= -146  & yi > -147; 
    end
    sumphv(ip).phv = zeros;
    sumphv(ip).phvstd = zeros;
end






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
        for i = 1:Nx
            for j=1:Ny
               GVv(i,j) = GVnow(i,j,ie,ip);
               GVv_std(i,j) = nanstd(GVnow(i,j,:,ip));
            end
        end
        GVv 
        if isSTD_weighted 
            GV(ip).wighted_region = (GVv_std <= 0.01) .* GV(ip).wighted_region;
        end
        GV_region = GVv .* GV(ip).wighted_region;   	 	
        GV_region
        ind = find(~isnan(GV_region) & GV_region ~= 0)
        evtavgphv(ip).phv(ie) = sum(GV_region(ind)) / size(ind,1);

        plot(periods(ip),evtavgphv(ip).phv(ie),'.b','markersize',15);
        hold on
        if (SW_type == 'L')
            load pa5pvL
        elseif ( SW_type == 'R')
            load pa5pvR    
        end
        plot(2*pi./w, phv,'r')

    end

    % plotting window setting
    xlim([10 110]);
    if (SW_type == 'L') 
        ylim([4.0 5.0]);
    elseif ( SW_type == 'R')
        ylim([4.0 4.2]);
    end
    xlabel('Period (s)');
    ylabel('Average Phase Velocity');
    gcarc(ie) = distance(eventinfo(ie).evla, eventinfo(ie).evlo, 9, -146);
    title(sprintf('Event %s dist %s',num2str(eventinfo(ie).id), num2str(gcarc(ie))));
    drawnow;


    if isfiguretoPS
       eventnamePS = [dispersion_curve_output_path,'/',num2str(eventinfo(ie).id),'.dispersion_curve.ps'];
       %print( 'depsc', eventname); doesn't work 
       %print -dpsc2 eventnamePS    doesn't work
        print('-dpsc2',eventnamePS)
    end


end

for ip = 1:length(periods)
    sumphv(ip).phv = nanmean(evtavgphv(ip).phv);
    ((evtavgphv(ip).phv-sumphv(ip).phv)/sumphv(ip).phv)*100
    sumphv(ip).phvstd = nanstd(evtavgphv(ip).phv);
end


figure(63)
clf


plot(periods,[sumphv.phv],'xb','linewidth',2);
hold on
errorbar(periods,[sumphv.phv],[sumphv.phvstd],'xb', 'linewidth',1);

%load allgrid
%plot(periods,phvallgrid,'xg')

if (SW_type == 'R')
    load pa5pvR
    plot(2*pi./w, phv,'r')
    hold on
	load age_52_110Myr_iso_pvR
	plot(period,phv,'k')

	hold on;
	load age_52_110Myr_aniso_pvR
	plot(period,phv,'k:')

elseif ( SW_type == 'L')
    load pa5pvL
    plot(2*pi./w, phv,'r')
    hold on
    load pa5pv1st
    plot(2*pi./w, phv,'r')
    hold on
end

%load age_52_110Myr_iso_DC
%plot(period,phv,'b')


xlim([10 110]);
if (SW_type == 'L')
    ylim([4.0 5.0]);
elseif ( SW_type == 'R')
    ylim([4.0 4.2]);
end

ylabel('Average Velocity (km/s)');
xlabel('Period (s)')
if isfiguretoPS
    AVG_dispersion_curvePS = [dispersion_curve_output_path,'/AVG.dispersion_curve.ps']
    %print -dpsc2 'AVG_dispersion_curve.ps' work
    print('-dpsc2',AVG_dispersion_curvePS)
end






%===================================
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
r = 0.01;
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(periods)
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    %h1=surfacem(xi,yi,avgphv(ip).GV_cor);
    h1=surfacem(xi,yi,avgphv(ip).GV.*GV(ip).wighted_region);
    % set(h1,'facecolor','interp');
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    avgv = nanmean(avgphv(ip).GV(:));
    if isnan(avgv)
        continue;
    end
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    load seiscmap
    colormap(seiscmap)
end
drawnow;

%==================================================

