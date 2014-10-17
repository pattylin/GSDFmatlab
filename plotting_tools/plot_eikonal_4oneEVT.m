function plot_eikonal_4oneEVT(EVT)
% Get information concnering station locations
%mk_stalist
load stalst;

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

load(['./eikonal/',num2str(EVT),'_eikonal_',comp,'.mat'])
load(['./eventmat/',num2str(EVT),'_',comp,'.mat'])


%CSnum = length(eventcs.CS); % CSnum = all pairs number for this events

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
    rays_temp = eventphv(ip).rays;
    w_temp = eventphv(ip).w;
   
    %temparray = vertcat(eventcs.CS.isgood);
   %goodind = find(temparray(:,ip) == 1 );
     %%
    h1=surfacem(xi,yi,eventphv(ip).GV);
    % set(h1,'facecolor','interp');
    %           load pngcoastline
    %           geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
    title(['Periods: ',num2str(periods(ip))],'fontsize',15)
    avgv = nanmean(eventphv(ip).GV(:));
    if isnan(avgv)
        continue;
    end
    r = 0.03;
    caxis([avgv*(1-r) avgv*(1+r)])
    colorbar
    load seiscmap
    colormap(seiscmap)
    for ir = 1:length(rays_temp)

        rayx = [rays_temp(ir,1);rays_temp(ir,3)];
        rayy = [rays_temp(ir,2);rays_temp(ir,4)];
        isgoodray = w_temp(ir);
        if isgoodray > 0
            plotm(rayx,rayy,'-k')
        else
            %plotm(rayx,rayy,'-','color', [0.5 0.5 0.5])
        end    
        plotm(xi,yi,'dk')
    end
     %title([num2str(length(goodind)) , 'Npairs in period  ',num2str(periods(ip))])


     

end
%  
% figure(46)
% clf
% N=3; M = floor(length(periods)/N)+1;
% clf
% for ip = 1:length(periods)
%     subplot(M,N,ip)
%     ax = worldmap(lalim, lolim);
%     set(ax, 'Visible', 'off')
%     if exist('drawlocal.m','file')
%          drawlocal
%     end
%     
%     
%     dt_refsta_temp = eventphv(ip).traveltime;
%     
%    
%     %temparray = vertcat(eventcs.CS.isgood);
%    %goodind = find(temparray(:,ip) == 1 );
%      %%
% 
%     for ir = 1:length(rays_temp)
%         
%         rayx = [rays_temp(ir,1);rays_temp(ir,3)];
%         rayy = [rays_temp(ir,2);rays_temp(ir,4)];
%         isgoodray = w_temp(ir);
%         if isgoodray > 0
%             plotm(rayx,rayy,'-k')
%         else
%             %plotm(rayx,rayy,'-','color', [0.5 0.5 0.5])
%         end    
%         %CS = eventcs.CS(ig);
%         %CSplot = eventcs.CSplot(ig);
%         %sta1 = CS.sta1;
%         %sta2 = CS.sta2;
%         %sta1_stnm = event.stadata(sta1).stnm;
%         %sta2_stnm = event.stadata(sta2).stnm;
% 
%         
%         %rayx = [event.stadata(sta1).stla;event.stadata(sta2).stla];
%         %rayy = [event.stadata(sta1).stlo;event.stadata(sta2).stlo];
%         %plotm(stla,stlo,'ok','markerfacecolor',GREEN);
%         %plotm(rayx,rayy,'-k');
%         plotm(xi,yi,'dk')
%      end
%      %title([num2str(length(goodind)) , 'Npairs in period  ',num2str(periods(ip))])
% 
% end
% 
% 
%  drawnow;

