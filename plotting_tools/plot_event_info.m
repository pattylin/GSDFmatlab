% For each event used plot the amplitude map, propagation anomaly map,
% and phase velocity so that we can figure out WHERE
% EVERYTHING IS GOING WRONG

% Will need to grab propagatoin anomaly and ray paths from eikonal and the
% rest from helmholtz arrays

% Nat-attack 2013

clear

% COLORS
BLUE=[44/255 77/255 143/255];
RED=[231/255 47/255 39/255];
GREEN=[19/255 166/255 50/255];

setup_parameters;
comp = parameters.component;
periods = parameters.periods;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize = parameters.gridsize;

% Set up the grid
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);
r = 0.1;

eikpath = './eikonal';
helmpath = './helmholtz';

load All_Stations.mat
load seiscmap.mat
% Get all the events used in the final stack from ueventlist.txt

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


% Start ploting
for iue = 30:length(uevent)
    phvmatfiles(iue) = dir([eikpath,'/',char(uevent(iue)),'_eikonal_',comp,'.mat']);
    temp = load([eikpath,'/',phvmatfiles(iue).name]);
    eventphv = temp.eventphv;
    disp(sprintf('Getting Information for %s',eventphv(1).id));
    stafile = sprintf('%12s_used_sta',char(uevent(iue)));
    
    % Convert things into dates so that only stations operating will be
    % plotted
    temp = sscanf(eventphv(1).id,'%4f %2f %2f %2f %2f');
    temp(6) = 0;
    ev_date = datenum(temp(1),temp(2),temp(3),temp(4),temp(5),temp(6));
    eg_stdate = datenum([2001 10 22 00 00 00]);
    ek_endate = datenum([2002 05 25 00 00 00]);
    clear stlas stlos
    if ev_date > eg_stdate && ev_date < ek_endate
        % Assign station locations to both EAGLE and EKBSE
        disp(sprintf('\t \t EAGLE and EKBSE active'));
        stlas = [stlas_ph1;stlas_EKSBE];
        stlos = [stlos_ph1;stlos_EKSBE];
        stnms = [stnms_ph1;stnms_EKSBE];
    elseif ev_date > ek_endate
        % Assign station locations to only EAGLE
        disp(sprintf('\t \t Only EAGLE active'));
        stlas = stlas_ph1;
        stlos = stlos_ph1;
        stnms = stnms_ph1;
    else
        % Assign station locations to only EKBSE
        disp(sprintf('\t \t Only EKBSE active'));
        stlas = stlas_EKSBE;
        stlos = stlos_EKSBE;
        stnms = stnms_EKSBE;
    end
    
    clear tempstla tempstlo
    for ista = 1:length(stlas)
        tempstla(ista) = str2num(sprintf('%6.4f',stlas(ista)));
        tempstlo(ista) = str2num(sprintf('%6.4f',stlos(ista)));
    end
    
    % First plot the ray paths ontop of the dynamic velocity
    disp(sprintf('\t Plotting Ray Paths'));
    N=3; M = floor(length(periods)/N) +1;
    clear pos
    f20=figure(20);
    clf
    fid2 = fopen(stafile,'w');
    if (fid2 == -1)
        error (['    Cannot open file: ', uevents]);
    end
    for ip = 1:length(periods)
        clear rays_temp
%         disp(ip)
        rays_temp = eventphv(ip).rays;
        subplot(M,N,ip)
        ax = worldmap(lalim,lolim);
        set(ax,'Visible','off')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[1 1 1],'EdgeColor','k')
        h1=surfacem(xi,yi,eventphv(ip).GV);
        colormap(seiscmap);
        avgv = nanmean(eventphv(ip).GV(:));
        if isnan(avgv)
            continue
        end
        caxis([avgv*(1-r) avgv*(1+r)])
        plotm(stlas,stlos,'ok','markerfacecolor',GREEN);
        %         hold on
        % set(h1,'facecolor','interp');
        % Open file that station names will be written to
        % Plot the individual rays and mark stations that are used red
        for ir = 1:length(rays_temp)
            rayx = [rays_temp(ir,1);rays_temp(ir,3)];
            rayy = [rays_temp(ir,2);rays_temp(ir,4)];
            plotm(rayx,rayy,'-k')
            plotm(rayx,rayy,'ok','markerfacecolor',RED);
            if ip == 4
                if ir == 1
                    fprintf(fid2,'Stations used in Cross Correlation for %12s\n'...
                        ,char(uevent(iue)));
                    disp(sprintf('\t Writing station names to file'));
                end
                rayx1 = str2num(sprintf('%6.4f',rayx(1)));
                rayx2 = str2num(sprintf('%6.4f',rayx(2)));
                ind1 = find(tempstla == rayx1);
                ind2 = find(tempstla == rayx2);
                rayn1 = stnms(ind1);
                rayn2 = stnms(ind2);
                fprintf(fid2,'%4s \t %4s\n',char(rayn1),char(rayn2));
            end
        end
        
    end
    fclose(fid2);
%     xxx
    mtit(sprintf('Ray Paths for %s',eventphv(1).id));
    drawnow;
    pos = get(f20,'Position');
    set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 660]);
    
    FIG1_id = [eventphv(1).id,'_rays.ps'];
    print(f20,'-dpsc',FIG1_id);
    %     xxx
    
    % Now plot propagation anomaly
    disp(sprintf('\t Plotting the propagation anomaly'));
    clear f21
    f21 = figure(21);
    clf
    [dtemp azi] = distance(xi(62),yi(62),eventphv(1).evla,eventphv(1).evlo);
    azi = azi - 180;
    azix = cosd(azi);
    aziy = sind(azi);
    for ip = 1:length(periods)
        subaxis(M,N,ip,'PR',0.01,'PL',0.01,'PT',0.02,'PB',0.02,'Margin',0.05);
        ax = worldmap(lalim,lolim);
        set(ax,'Visible','off')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[1 1 1],'EdgeColor','k')
        surfacem(xi,yi,eventphv(ip).azi_diff)
        hold on
        quiverm(xi(62),yi(62),azix,aziy,'k');
        colormap(seiscmap);
        cbar_axis = colorbar();
        set(get(cbar_axis,'xlabel'),'String', 'degree');
    end
    pos = get(f21,'Position');
    set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 660]);
    mtit(sprintf('Propagation Direction Anomaly for %s',eventphv(1).id));
    drawnow;
    FIG2_id = [eventphv(1).id,'_propanom.ps'];
    print(f21,'-dpsc',FIG2_id);
    % %     xxx
    
    %% Now grab information from helmholtz arrays
    helmmatfiles(iue) = dir([helmpath,'/',char(uevent(iue)),'_helmholtz_',comp,'.mat']);
    temp = load([helmpath,'/',helmmatfiles(iue).name]);
    helmholtz = temp.helmholtz;
    %     disp(helmholtz(1).id);
    
    % Plot structural phase velocity
    disp(sprintf('\t Plotting Structural Phase Velocity'));
    f22 = figure(22);
    clf
    for ip = 1:length(periods)
        %         subplot(M,N,ip)
        subaxis(M,N,ip,'PR',0.01,'PL',0.01,'PT',0.02,'PB',0.02,'Margin',0.05);
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[1 1 1],'EdgeColor','k')
        h1=surfacem(xi,yi,helmholtz(ip).GV_cor);
        %         set(h1,'facecolor','interp');
        avgv = nanmean(helmholtz(ip).GV(:));
        colorbar
        load seiscmap
        colormap(seiscmap)
        if isnan(avgv)
            caxis([3 4])
        else
            caxis([avgv*(1-r) avgv*(1+r)])
        end
    end
    pos = get(f22,'Position');
    set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 660]);
    mtit(sprintf('Structural Phase Velocity for %s',helmholtz(1).id));
    drawnow;
    FIG3_id = [eventphv(1).id,'_phasev.ps'];
    print(f22,'-dpsc',FIG3_id);
    
    % Plot amplitude map
    disp(sprintf('\t Plotting the Amplitude Map'));
    f23 = figure(23);
    clf
    for ip = 1:length(periods)
        %         subplot(M,N,ip)
        subaxis(M,N,ip,'PR',0.01,'PL',0.01,'Margin',0.05);
        ax = worldmap(lalim,lolim);
        set(ax,'Visible','off')
        h1=surfacem(xi,yi,helmholtz(ip).ampmap);
        hold on
        %         plotm(stlas_EKSBE,stlos_EKSBE,'ok','markerfacecolor',GREEN);
        hold on
        colorbar
        colormap(seiscmap)
    end
    pos = get(f23,'Position');
    set(gcf, 'PaperOrientation', 'landscape','PaperPositionMode','auto',...
        'Position',[pos(1) pos(2) 700 660]);
    mtit(sprintf('Amplitude Mape for %s',helmholtz(1).id));
    drawnow;
    FIG4_id = [eventphv(1).id,'_amplitude.ps'];
    print(f23,'-dpsc',FIG4_id);
%     pause;
end


