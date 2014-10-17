% Plot histogram of observations at each period
% Reads in information from ueventlist.txt that is created by stack_helm.m
clear

setup_parameters;

periods = parameters.periods;

uevents = 'ueventlist.txt';
fid = fopen(uevents,'r');
if (fid == -1)
    error (['    Cannot open file: ', uevents]);
end
f = textscan(fid,'%f %s');
fclose(fid);

per_obs = f{1};

figure(61)
clf
hist(per_obs,periods);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')
xlabel('Period (s)');
title('Number of Observations per Period');