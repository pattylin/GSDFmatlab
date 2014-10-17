% Plot phase velocity at a given location against the azimuth

clear
chooselocation = 0;
RED = [231/255 47/255 39/255];
GREEN = [19/255 166/255 50/255];
BLUE = [44/255 77/255 143/255];
PEACH = [245/255 223/255 181/255];

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
r = 0.10;

% Needed to plot the standard deviation
load(['helmholtz_stack_',comp,'.mat']);

Nx = length(avgphv(1).xnode);
Ny = length(avgphv(1).ynode);
xi = avgphv(1).xi;
yi = avgphv(1).yi;

helmpath = './helmholtz';

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


%% Plot at multiple locations
% It could also be done such that you create a vector of locations you want
% to look at and then loop through them only rather than the entire grid

disp('Let us start walking through the array');
for ix = 1:(Nx*Ny)
    
    uxi = ix;
    
    % Skip locations along the limits of the grid
    if xi(uxi) == lalim(1) || xi(uxi) == lalim(2) || xi(uxi) > 11
        continue
    elseif yi(uxi) == lolim(1) || yi(uxi) == lolim(2)
        continue
    end
    disp(sprintf('Location %s N %s E',num2str(xi(uxi)),num2str(yi(uxi))));
    boxlat = [(xi(uxi)-.25) (xi(uxi)+.25) (xi(uxi)+.25) (xi(uxi)-.25) (xi(uxi)-.25)];
    boxlon = [(yi(uxi)-.25) (yi(uxi)-.25) (yi(uxi)+.25) (yi(uxi)+.25) (yi(uxi)-.25)];
    
    
    figure(42)
    clf
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 .5 .5])
    h1=surfacem(xi,yi,avgphv(4).GV_cor);
    hold on
    %     plotm(xi(uxi),yi(uxi),'v')
    box = [boxlat' boxlon'];
    plotm(box,'-w','linewidth',2);
    hold on
    plotm(box,'.k','markersize',21);
    hold off;
    colorbar
    load seiscmap
    colormap(seiscmap)
    
    % Set up all the pretty colors
    CC = hsv(length(periods));
    
    %% Get information concerning the phase velocity and azimuth
    clear h k
    figure(43)
    clf
    if ~isnan(avgphv(5).GV_cor_std(uxi)) && ~isnan(avgphv(5).GV_cor(uxi))
        tempazi = [-180;180];
        GVstd_p = avgphv(5).GV_cor_std(uxi) + avgphv(5).GV_cor(uxi);
        GVstd_n = -avgphv(5).GV_cor_std(uxi) + avgphv(5).GV_cor(uxi);
        std_temp = [GVstd_p;GVstd_p];
        k = area(tempazi,std_temp);
        set(k,'BaseValue',GVstd_n,'FaceColor',PEACH,...
            'LineStyle',':');
        
        hold on
    end
    
    clear evazi evphv
    for ip=1:length(periods)
       
        for iue = 1:length(uevent)
            clear dtemp azi
            helmmatfiles(iue) = dir([helmpath,'/',char(uevent(iue)),'_helmholtz_',comp,'.mat']);
            temp = load([helmpath,'/',helmmatfiles(iue).name]);
            helmholtz = temp.helmholtz;
            if isnan(helmholtz(ip).GV_cor(uxi))==1;
                %                 evazi(ip).event(iue) = NaN;
                %                 evphv(ip).event(iue) = NaN;
                continue
            else
                %         disp(helmholtz(1).id);
                [dtemp azi] = distance(xi(62),yi(62),helmholtz(1).evla,helmholtz(1).evlo);
                evazi(ip).event(iue) = azi - 180;
                evphv(ip).event(iue) = helmholtz(ip).GV_cor_std(uxi);
            end
        end % end of events loop
        
        
        if exist('evazi')==0 || length(evphv) < ip
            continue
        else
            temp = find(evphv(ip).event == 0);
            evphv(ip).event(temp) = NaN;
            evazi(ip).event(temp) = NaN;
            clear temp
            temp = find(evphv(ip).event > 0);
            disp(sprintf('%s Events Used Here with %s non-NaN results for %s s',...
                num2str(length(evphv(ip).event)),num2str(length(temp)),num2str(periods(ip))));
            
            
            % Plot the average structural velocity dispersion for the chosen location
            h(ip) = plot(evazi(ip).event(~isnan(evazi(ip).event))...
                ,evphv(ip).event(~isnan(evphv(ip).event)),'.r','markersize',...
                20,'color',CC(ip,:));
            xlabel('Azimuth (degree)')
            ylabel('Velocity (km/s)');
            xlim([-180 180]);
            %             ylim([3 4.5]);
            title(sprintf('Structural Phase Velocity at %s N %s E',num2str(xi(uxi)),...
                num2str(yi(uxi))));
            hold on
            
        end
        
        hold on;
    end % end of periods loop
    hold off;
    if exist('evazi')==0
        continue
    else
        htemp = find(h ~= 0);
        if length(htemp) == 8
            legend(h,'20','25','32','40','50','60','80','100');
        end
        drawnow;
%         pause
    end
end
