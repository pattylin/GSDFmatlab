clear all;

setup_parameters;

periods = parameters.periods;
comp = parameters.component;


CSpath = './CSmeasure/';
OPWpath = './OPW/';
eventcsfiles = dir([CSpath,'/*_',comp,'.mat']);
num_evt = length(eventcsfiles);
colorsR = linspace(1,0,num_evt) ;
colorsG = linspace(0,1,num_evt) ;
colorsB = linspace(0,1,num_evt) ;
colors =  [colorsR ; colorsG ; colorsB] ;

% read in bad event list, if existed
if exist('badevt.lst')
    badevts = textread('badevt.lst','%s');
    disp('Found Bad event')
    disp(badevts)
end

figure(101)
clf
load([CSpath,eventcsfiles(1).name])

refsta = 'B13'
indrefsta = find(strncmp(refsta,eventcs.stnms,3));
for ip = 1:length(periods)


    phv_array = [];
    theta_array = [];
	for ie = 1:length(eventcsfiles)
		%if ~exist('badevt.lst') ||
		%isempty(find(strcmp(eventcsfiles(ie).name(1:12),badevts), 1)) %
		%works too. 
        clear sta_phv_evt
        if ~exist('badevt.lst') || isempty(find(strncmpi(eventcsfiles(ie).name,badevts,12), 1))
            %eventcsfiles(ie).name;
            %[sta_phv sta_theta sta_gcarc] = OPW_inv_4oneEVT([CSpath,eventcsfiles(ie).name],ip); 
       
            matfilename = [OPWpath,'/',eventcsfiles(ie).name(1:12),'_OPW_',comp,'.mat'];
            load(matfilename);
            sta_phv = ArraySurface(ip).sta_phv(indrefsta);
            sta_azi = ArraySurface(ip).sta_theta(indrefsta);
            sta_gcarc = ArraySurface(ip).sta_gcarc(indrefsta);
        
            ind = find(abs(sta_phv-nanmean(sta_phv))/nanmean(sta_phv)*100 >=5 );
            if ~isempty(ind)
                sta_phv(ind) = [];
                sta_azi(ind)= [];
            end
        end        
        
        phv_array =[phv_array,sta_phv];
        theta_array = [theta_array,sta_azi];
		%phv_array =[phv_array,nanmean(sta_phv)];
        %theta_array = [theta_array, nanmean(sta_azi)]; 
    end
    
    ind = find(~ isnan(phv_array));

    theta_array = theta_array - 180.00;
    [para fiterr]=fit_azi_anisotropy(theta_array(ind),phv_array(ind));
    parastd=confint(para);
    number_phv = length(phv_array);
    
    isophv_array=para.a;
    isophv_array_std=parastd(2,1)-parastd(1,1);
    aniso_array_strength=para.d;
    aniso_array_azi=para.e;

    
    
    if para.e > 180
        aniso_array_azi=para.e-180;
    elseif para.e < 0
        aniso_array_azi=para.e+180;
    end
    
    
    aniso_strength_std =parastd(2,2)-parastd(1,2);
    aniso_azi_std      =parastd(2,3)-parastd(1,3);
    
    
    OPWSurface(ip).isophv = isophv_array;
    OPWSurface(ip).isophv_std = isophv_array_std;
    OPWSurface(ip).aniso_strength = aniso_array_strength;
    OPWSurface(ip).aniso_azi = aniso_array_azi;
    OPWSurface(ip).aniso_strength_std = aniso_strength_std;
    OPWSurface(ip).aniso_azi_std = aniso_azi_std;
    OPWSurface(ip).parameters = parameters;
    OPWSurface(ip).numer_phvmeasurment = length(ind);
    
    plot_position = [1 3 5 7 2 4 6 8 ];
    
    subplot(4,2,plot_position(ip))
    hold on
    allazi = -200:200;
    plot(allazi,(para.a*(1+para.d*cosd(2*(allazi-para.e)))-mean(phv_array(ind)))/mean(phv_array(ind))*100,'r','linewidth',2)


    plot(theta_array(ind),(phv_array(ind)-mean(phv_array(ind)))/mean(phv_array(ind))*100,'x','linewidth',2); hold on;

    %length(azi)
    %max(abs((sta_phv_array-mean(sta_phv_array))/mean(sta_phv_array)*100))
    %[azi (phV-mean(phV))/mean(phV)*100]
    xlim([-180 180])
    set(gca,'xtick',[-180 -135 -90 -45 0 45 90 135 180])
    title(['Period= ',num2str(periods(ip)),' (s)'])
    ylim([-4 4])
    set(gca,'ytick',[-4:2:4])
    %pause;
    hold on;

    


   
end

for ip = 1:length(periods)
    phv_1(ip) = OPWSurface(ip).isophv;
    phv_1_std(ip) = OPWSurface(ip).isophv_std/sqrt(OPWSurface(ip).numer_phvmeasurment);
    
    aniso_strength_1(ip) = OPWSurface(ip).aniso_strength;
    
    aniso_strength_1_std(ip) = OPWSurface(ip).aniso_strength_std/sqrt(OPWSurface(ip).numer_phvmeasurment);
    
    aniso_azi_1(ip) = OPWSurface(ip).aniso_azi;
    
    aniso_azi_1_std(ip) = OPWSurface(ip).aniso_azi_std/sqrt(OPWSurface(ip).numer_phvmeasurment);
    
end

figure(201)
clf
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
xlim([10 110]);
ylim([4.0 4.2]);
legend('isophv','pa5', 'age52to110iso','age52to110aniso')
ylabel('Phase Velocity (km/s)');



subplot(3,2,3)
plot(periods,aniso_strength_1*100,'xb','linewidth',2)
errorbar(periods,aniso_strength_1*100,aniso_strength_1_std*100,'xb', 'linewidth',1);
ylabel('amplitude');
xlim([10 110]);
ylim([0 3]);



subplot(3,2,5)
plot(periods,aniso_azi_1,'xb','linewidth',2)
errorbar(periods,aniso_azi_1,aniso_azi_1_std,'xb', 'linewidth',1);
ylabel('fast direction');
xlim([10 110]);
ylim([0 180]);
set(gca,'ytick',[0 45 90 135 180])
xlabel('Period (s)')


