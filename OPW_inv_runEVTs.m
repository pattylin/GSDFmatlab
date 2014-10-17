clear all;

setup_parameters;

periods = parameters.periods;
comp = parameters.component;


CSpath = './CSmeasure/';
eikonalpath = './eikonal/';
OPW_outputpath = './OPW/';

if ~exist(OPW_outputpath)
	mkdir(OPW_outputpath);
end

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

iscal = 1;
if iscal 
    for ie = 1: length(eventcsfiles)

        eventcsfiles(ie).name
        %if ~exist('badevt.lst') ||
            %isempty(find(strcmp(eventcsfiles(ie).name(1:12),badevts), 1)) %
            %works too. 
        if ~exist('badevt.lst') || isempty(find(strncmpi(eventcsfiles(ie).name,badevts,12), 1))
            for ip = 1:length(periods)


                [sta_phv sta_theta sta_gcarc] = OPW_inv_4oneEVT([eventcsfiles(ie).name(1:12)],ip);

                ArraySurface(ip).sta_phv = sta_phv;
                ArraySurface(ip).sta_theta = sta_theta;
                ArraySurface(ip).sta_gcarc = sta_gcarc; 
                

            end


        end

        matfilename = [OPW_outputpath,'/',eventcsfiles(ie).name(1:12),'_OPW_',comp,'.mat'];
        save(matfilename,'ArraySurface');
        disp(['Save the result to: ',matfilename])

    end
end



isplot = 1;

figure(101)
clf;
if isplot 
    for ip = 1:length(periods)


        for ie = 1: length(eventcsfiles)
            if ~exist('badevt.lst') || isempty(find(strncmpi(eventcsfiles(ie).name,badevts,12), 1))
                matfilename = [OPW_outputpath,'/',eventcsfiles(ie).name(1:12),'_OPW_',comp,'.mat'];
                load(matfilename);
                size(ArraySurface(ip).sta_phv);
                sta_phv = ArraySurface(ip).sta_phv(9);
                sta_azi = ArraySurface(ip).sta_theta(9);
                sta_gcarc = ArraySurface(ip).sta_gcarc(9);
                clear ind;
                
                ind = find(abs(sta_phv-nanmean(sta_phv))/nanmean(sta_phv)*100 >=5 );
                if ~isempty(ind)
                    sta_phv(ind) = [];
                    sta_gcarc(ind)= [];
                end
                %ind
                
                
                %sta_phv(~isnan(sta_phv));

                 plot_position = [1 3 5 7 2 4 6 8 ];
    
                subplot(4,2,plot_position(ip))

                plot(sta_azi,sta_phv,'.','color',colors(:,ie));hold on; 
                



            end
        end

        %xlabel('gcarc (degree)');
    end
end

 %plotting
isplot = 1;
if isplot 
    figure(102)
    clf
    for ip = 1:length(periods)

        plot_position = [1 3 5 7 2 4 6 8 ];

        subplot(4,2,plot_position(ip))

        for ie = 1: length(eventcsfiles)
            if ~exist('badevt.lst') || isempty(find(strncmpi(eventcsfiles(ie).name,badevts,12), 1))
                matfilename = [OPW_outputpath,'/',eventcsfiles(ie).name(1:12),'_OPW_',comp,'.mat'];
                load(matfilename);
                sta_phv = ArraySurface(ip).sta_phv(9);
                sta_azi = ArraySurface(ip).sta_theta(9);
                sta_gcarc = ArraySurface(ip).sta_gcarc(9);
                clear ind;
                
                ind = find(abs(sta_phv-nanmean(sta_phv))/nanmean(sta_phv)*100 >=5 );
                if ~isempty(ind)
                    sta_phv(ind) = [];
                    sta_gcarc(ind)= [];
                end
                
                ind = find(~isnan(sta_phv));
                
                

                
               
                plot(sta_gcarc(ind),sta_phv(ind),'.','color',colors(:,ie));hold on; 
                ylabel('phv (grid-search slowness vector)');hold on;
                xlim([30 100]); hold on;

                %subplot(2,1,2)
                %plot(sta_gcarc,(sta_phv-nanmean(sta_phv))*100,'.','color',colors(:,ie)); hold on;
                %%ylim([-15 15]); hold on;

            end
        end

        %xlabel('gcarc (degree)');
    end
end
