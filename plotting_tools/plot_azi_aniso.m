addpath('/Users/pappi/Research/PGS/MATLABdatalib')

% (4,4) = grid(25)%(2,2)

periods = parameters.periods;
for ip=1:length(periods)

    phv_1(ip) = avgphv_aniso(ip).isophv(2,2);
    %phv_1_std(ip) = avgphv_aniso(ip).isophv_std(2,2);  
    phv_1_std(ip) = avgphv_aniso(ip).isophv_std(2,2)/sqrt(avgphv_aniso(ip).numer_phvmeasurment(2,2));  
    aniso_strength_1(ip) = avgphv_aniso(ip).aniso_strength(2,2);
    %aniso_strength_1_std(ip) = avgphv_aniso(ip).aniso_strength_std(2,2);
    aniso_strength_1_std(ip) = avgphv_aniso(ip).aniso_strength_std(2,2)/sqrt(avgphv_aniso(ip).numer_phvmeasurment(2,2));
    aniso_azi_1(ip) = avgphv_aniso(ip).aniso_azi(2,2);
    %aniso_azi_1_std(ip) = avgphv_aniso(ip).aniso_azi_std(2,2);
    aniso_azi_1_std(ip) = avgphv_aniso(ip).aniso_azi_std(2,2)/sqrt(avgphv_aniso(ip).numer_phvmeasurment(2,2));
end



figure(101)
%clf
hold on;
subplot(3,2,1)
plot(periods,phv_1,'xb','linewidth',2)
errorbar(periods,phv_1,phv_1_std,'xb', 'linewidth',1);
hold on
load pa5pvR

plot(2*pi./w, phv,'r')
hold on;
load age_52_110Myr_iso_pvR
plot(period,phv,'k')

hold on;
load age_52_110Myr_aniso_pvR
plot(period,phv,'k:')
%xlim([10 110]);
xlim([10 180])
ylim([4.0 4.3]);
legend('isophv','pa5', 'age52to110iso','age52to110aniso')
ylabel('Phase Velocity (km/s)');



subplot(3,2,3)
hold on;
peak2peak_amp = 2 .*aniso_strength_1
plot(periods, peak2peak_amp*100,'xb','linewidth',1)
errorbar(periods,peak2peak_amp*100,aniso_strength_1_std*100,'xb', 'linewidth',1);
hold on;
ylabel('peak-to-peak amp (%)');
%xlim([10 110]);
xlim([10 180]);
ylim([0 6]);



subplot(3,2,5)
hold on;
plot(periods,aniso_azi_1,'xb','linewidth',1); 
errorbar(periods,aniso_azi_1,aniso_azi_1_std,'xb', 'linewidth',1);
hold on;
ylabel('fast direction');
xlim([10 180]);
ylim([0 180]);
set(gca,'ytick',[0 45 90 135 180])
xlabel('Period (s)')





psfile = ['azimutah_aniso.ps']
print('-dpsc2',psfile)
