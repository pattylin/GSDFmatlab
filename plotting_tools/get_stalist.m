% Script to output used station pairs
% NJA 2013

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

CSpath = './CSmeasure';

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


for iue = 1:length(uevent)
    CSmat(iue) = dir([CSpath,'/',char(uevent(iue)),'_CS_',comp,'.mat']);
    temp = load([CSpath,'/',CSmat(iue).name]);
    eventcs = temp.eventcs;
    stnms = eventcs.stnms;
    stafile = sprintf('%12s_used_sta',char(uevent(iue)));
    fid2 = fopen(stafile,'w');
    if (fid2 == -1)
        error (['    Cannot open file: ', uevents]);
    end
    
    for ics = 1:length(eventcs.CS)
        CS = eventcs.CS(ics);
        good = 0;
        for ip = 1:length(periods)
            if good == 1;
                continue
            elseif CS.isgood(ip) == 1 && good == 0;
                stn1 = stnms(CS.sta1);
                stn2 = stnms(CS.sta2);
                fprintf(fid2,'%4s  %4s\n',char(stn1),char(stn2));
                good = 1;
            else
                continue
            end
        end
        
        
    end
end
fclose(fid2);

                
    
    
