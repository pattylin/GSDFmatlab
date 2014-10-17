% Program to plot the best alpha distribution for all frequency bands.

clear;

setup_parameters;
phase_v_path = './helmholtz/'
comp = parameters.component;
periods = parameters.periods;

phvmatfiles = dir([phase_v_path,'/*_helmholtz_',comp,'.mat']);

for ie = 1:length(phvmatfiles)
    disp(phvmatfiles(ie).name);
    temp = load([phase_v_path,phvmatfiles(ie).name]);
    helmholtz = temp.helmholtz;
    alpha(:,ie) = [helmholtz.bestalpha]';
    
    for ip=1:length(periods)
        dist(ip,ie) = [helmholtz(ip).dist]';
        dist(ip,ie) = dist(ip,ie);
        %phv(ip,ie) = [helmholtz(ip).GV_cor(68)]';
        phv(ip,ie) = helmholtz(ip).GV_cor(68);
    end
end

figure(12)
clf
N=3; M = floor(length(periods)/N)+1;
for ip=1:length(periods);
    subplot(M,N,ip)
    hist(alpha(ip,:));
    title(sprintf(' %s s',num2str(periods(ip))));
end


N=3; M = floor(length(periods)/N)+1;
figure(13)
clf
for ip=1:length(periods);
    phv_red(ip,:) = phv(ip,:)./nanmean(phv(ip,:));
    subplot(M,N,ip)
    plot(dist(ip,:),phv_red(ip,:),'kx','markersize',10,'linewidth',2)
    ylabel('Normalized Phase Velocity');
    ylim([0.85 1.25]);
    %ylabel('Phase Velocity (km/s)');
    xlabel('Epicentral Distance (deg)');
    title(sprintf(' %s s reduced by %s',num2str(periods(ip)),num2str(nanmean(phv(ip,:)))));
end
drawnow;
