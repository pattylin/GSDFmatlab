% Script to plot events used in GSDF
% Natalie Accardo
clear all

RED = [231/255 47/255 39/255];
GREEN = [19/255 166/255 50/255];
BLUE = [44/255 77/255 143/255];

%load All_Stations.mat

dbpath = './sacdata/';
eventfile = 'ueventlist.txt';
outpath = './eventmat/';

setup_parameters;
comp = parameters.component;

buff = 5; % Decides size of box to search in

uevents = 'ueventlist.txt';
fid = fopen(uevents,'r');
if (fid == -1)
    error (['    Cannot open file: ', uevents]);
end
f = textscan(fid,'%f %s');
fclose(fid);

temp=f{2};
count = 0;
for ie = 1:(length(temp)-1)
    if strcmp(temp(ie),temp(ie+1)) ==1
        continue
    else
        count = count+1;
        eventids(count) = temp(ie);
    end
end


for ie=1:length(eventids)
    matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
    %     disp(matfilename)
    load(matfilename)
    evid(ie,1:12) = event;
    evla(ie) = event.evla;
    evlo(ie) = event.evlo;
    evmag(ie) = event.mag;
end

figure(50)
clf
% set(gcf,'FontSize',16)
ax = worldmap('World');
setm(ax, 'Origin', [0 180 0])
p = findobj(ax,'type','patch');
land = shaperead('landareas','UseGeocoords',true);
geoshow(ax, land, 'FaceColor', [1 0.75 0.01])
set(p,'FaceColor',[.49 .73 .82]);
hold on
%plotm(stlas_EKSBE,stlos_EKSBE,'v','markerfacecolor',BLUE,'markersize',7,'markeredgecolor',BLUE);
%plotm(stlas_ph1,stlos_ph1,'v','markerfacecolor',BLUE,'markersize',7,'markeredgecolor',BLUE);
for ie=1:length(eventids)
    if evmag(ie) >= 6 && evmag(ie) < 7
        plotm(evla(ie),evlo(ie),'o','markerfacecolor',GREEN,'markersize',12,'markeredgecolor',GREEN);
        hold on
    elseif evmag(ie) >= 7 && evmag(ie) < 8
        plotm(evla(ie),evlo(ie),'o','markerfacecolor',RED,'markersize',15,'markeredgecolor',RED);
        hold on
    else
        plotm(evla(ie),evlo(ie),'bv','markerfacecolor','b');
    end
    hold on
end
title('Events Used','FontSize',16,'fontweight','bold');
hold off
ans = input('Do you want to select an Earthquake? Y/N [Y]:','s');
if strcmp(ans,'Y') ==1
    count = 0;
    while count < 10
        [x,y]=inputm(1);
        
        A = find(evla > (x-buff) & evla < (x+buff) & evlo > (y-buff) & evlo < (y+buff));
        disp(sprintf('Event %s Lat %3.2f Lon %4.2f Mag %3.2d',char(eventids(A)),evla(A),evlo(A),evmag(A)));
        ans2 = input('Do you want to find another EQ? Y/N [Y]:','s');
        if strcmp(ans2,'Y') == 1
            count = count +1;
            continue
        else
            count = 11;
            disp('This means that we are done for now ...');
        end
    end
else
    disp('See you next time then')
end
