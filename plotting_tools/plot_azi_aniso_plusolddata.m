addpath('/Users/pappi/Research/PGS/MATLABdatalib')

% (4,4) = grid(25)%(2,2)
parameters.periods = [20 25 32 40 50 60 80 100];
load ../GSDF_Rayleigh/gridsize3.0_maxdist600/eikonal_stack_aniso_LIZ.mat
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
plot(periods,phv_1,'xc','linewidth',2)
errorbar(periods,phv_1,phv_1_std,'xc', 'linewidth',1);
hold on

xlim([10 180])
ylim([4.0 4.3]);



subplot(3,2,3)
hold on;
peak2peak_amp = 2 .*aniso_strength_1
plot(periods, peak2peak_amp*100,'xc','linewidth',1)
errorbar(periods,peak2peak_amp*100,aniso_strength_1_std*100,'xc', 'linewidth',1);
hold on;
xlim([10 180]);
ylim([0 6]);



subplot(3,2,5)
hold on;
plot(periods,aniso_azi_1,'xc','linewidth',1);
errorbar(periods,aniso_azi_1,aniso_azi_1_std,'xc', 'linewidth',1);
hold on;

xlim([10 180]);
ylim([0 180]);





psfile = ['azimutah_aniso.ps']
print('-dpsc2',psfile)
