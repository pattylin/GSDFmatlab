% Generate mat file for station names and locations
% NJA, 2013 LDEO
clear;

isfigure = 1;
is_overwrite = 0;

eventmatpath = './eventmat/';
outfile = 'stalst.mat'

% Setup parameters
setup_parameters

% Setup Error Codes for Bad data
setup_ErrorCode

periods = parameters.periods;
comp = parameters.component;
count = 0;
matfiles = dir([eventmatpath,'/*_',comp,'.mat']);

for ie = 1:length(matfiles)
    
    clear event
    % read in the events information
    temp = load([eventmatpath,matfiles(ie).name]);
    event = temp.event;
    disp(event.id)
    for ista = 1:length(event.stadata)
        if exist('stla','var') == 0
            disp('Create the Array')
            count = count+1;
                stla(count) = event.stadata(ista).stla;
                stlo(count) = event.stadata(ista).stlo;
                stnm(count) = cellstr(event.stadata(ista).stnm);
            continue
        else
%            A = find(event.stadata(ista).stla == lat & ...
%                event.stadata(ista).stlo == lon);
            temp = strcmp(event.stadata(ista).stnm,stnm);
            A = find(temp == 1);
            if A > 0
                continue
            else
                count = count+1;
                stla(count) = event.stadata(ista).stla;
                stlo(count) = event.stadata(ista).stlo;
                stnm(count)= cellstr(event.stadata(ista).stnm);
            end
        end
    end
end

save(outfile,'stla','stlo','stnm');
