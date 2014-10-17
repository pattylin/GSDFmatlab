% Plot structural phase velocity and standard deviation across all periods
% at a given location
% by accardo

clear
chooselocation = 0;
RED = [231/255 47/255 39/255];
GREEN = [19/255 166/255 50/255];
BLUE = [44/255 77/255 143/255];

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
r = 0.10;

load helmholtz_stack_BHZ.mat
Nx = length(avgphv(1).xnode);
Ny = length(avgphv(1).ynode);
xi = avgphv(1).xi;
yi = avgphv(1).yi;

if chooselocation == 1
    % Choose the location you want to look at
    ulat = 8;
    ulon = -146;
    basem = zeros([Nx Ny]);
    
    % Plot a map showing the location you have chosen
    tempx = find(xi == ulat);
    tempy = find(yi == ulon);
    for ix = 1:length(tempx)
        for iy = 1:length(tempy)
            if tempx(ix) == tempy(iy)
                basem(tempx(ix)) = 1;
                uxi = tempx(ix);
                uyi = tempy(iy);
            end
        end
    end
end

% If you only want to look at one location turn off the for loop and next
% two lines

for ix = 1:(Nx*Ny)
    uxi = ix;
    disp(xi(uxi));
    disp(yi(uxi));
    
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
    
    % Plot the average structural velocity dispersion for the chosen location
    figure(43)
    clf
    for ip = 1:length(periods)
        if isnan(avgphv(ip).GV_cor(uxi))
            plot(periods(ip),3.5,'xb','markersize',10,'linewidth',2)
        end
        plot(periods(ip),avgphv(ip).GV_cor(uxi),'.r','markersize',10)
        errorbar(periods(ip),avgphv(ip).GV_cor(uxi),avgphv(ip).GV_cor_std(uxi))
        hold on
    end
    ylim([3 4])
    xlabel('Period (s)')
    ylabel('Velocity (km/s)');
    title(sprintf('Structural Phase Velocity at %s N %s E',num2str(xi(uxi)),...
        num2str(yi(uxi))));
    drawnow;
    
    % Plot standard deviation for the chosen location
%     figure(44)
%     clf
%     for ip = 1:length(periods)
%         plot(periods(ip),avgphv(ip).GV_cor_std(uxi),'.k','markersize',20,'color',GREEN);
%         hold on
%     end
%     drawnow;
%    pause
end
