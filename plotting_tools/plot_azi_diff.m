% Program to Plot the average propagation direction anomaly
% You'll want to run this with isssaved=0 only when the eikonal mat file
% has already been appended with the azimuth information

clear;

issaved=1;
phase_v_path = './eikonal/'
r = 0.10;

% COLORS
BLUE=[44/255 77/255 143/255];
RED=[231/255 47/255 39/255];
GREEN=[19/255 166/255 50/255];

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
mincsnum = parameters.mincsnum;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
min_event_num = parameters.min_event_num;
err_std_tol = parameters.err_std_tol;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

for ip=1:length(periods)
    avgphv(ip).azi_avg = zeros(Nx,Ny);
end

phvmatfiles = dir([phase_v_path,'/*_',comp,'.mat']);
if (issaved==1)
    for ip = 1:length(periods)
        disp(sprintf('Period : %s',num2str(periods(ip))))
        for in = 1:(Nx*Ny)
            clear temp_azi
            for ie = 1:length(phvmatfiles)
                temp = load([phase_v_path,phvmatfiles(ie).name]);
                eventphv = temp.eventphv;
                %             disp(eventphv(1).id);
                temp_azi(ie) = eventphv(ip).azi_diff(in);
            end
            avgphv(ip).azi_avg(in) = nanmean(temp_azi);
            
        end
        azi(ip).avg = avgphv(ip).azi_avg;
    end
    
    save(['eikonal_stack_',comp,'.mat'],'azi','-append');
else
    load eikonal_stack_BHZ.mat
end
N=3; M = floor(length(periods)/N)+1;

figure(100)
clf
load All_Stations.mat
for ip = 1:length(periods)
    subplot(M,N,ip)
    ax = worldmap(lalim,lolim);
    set(ax,'Visible','off');
%         h1=surfacem(xi,yi,avgphv(ip).azi_avg);
    h1=surfacem(xi,yi,azi(ip).avg);
    set(h1,'facecolor','interp')
    plotm(stlas_EKSBE,stlos_EKSBE,'ok','markerfacecolor',GREEN);
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    colorbar
    load seiscmap
    colormap(seiscmap)
end
P = mtit('Summed Propagation Direction Anomaly');
