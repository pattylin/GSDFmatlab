clear all;

setup_parameters;

periods = parameters.periods;
comp = parameters.component;
load stalst;

figure(101)
clf



CSpath = './CSmeasure/';

eventcsfiles = dir([CSpath,'/*_',comp,'.mat']);
num_evt = length(eventcsfiles);
colorsR = linspace(1,0,num_evt) ;
colorsG = linspace(0,1,num_evt) ;
colorsB = linspace(0,1,num_evt) ;
colors =  [colorsR ; colorsG ; colorsB] ;


for ip = 4:4
    
    figure(101)
    clf
    for ie = 1:length(eventcsfiles)
       clear eventcs
       load([CSpath,eventcsfiles(ie).name])
       
       avgphv_evt = eventcs.avgphv(ip);
       for ics = 1: length(eventcs.CS)

            ddist(ics) = eventcs.CS(ics).ddist;
            dtp(ics) = eventcs.CS(ics).dtp(ip);
            dtg(ics) = eventcs.CS(ics).dtg(ip);
            isgood(ics) = eventcs.CS(ics).isgood(ip);

            if isgood(ics) == 1
                sta1_stnm = stnm(eventcs.CS(ics).sta1);
                sta2_stnm = stnm(eventcs.CS(ics).sta2);
                sta1_stla = stla(eventcs.CS(ics).sta1);
                sta2_stla = stla(eventcs.CS(ics).sta2);
                sta1_stlo = stlo(eventcs.CS(ics).sta1);
                sta2_stlo = stlo(eventcs.CS(ics).sta2);
                arclen(ics) = distance( eventcs.evla, eventcs.evlo, (sta1_stla+ sta2_stla)/2, (sta1_stlo+ sta2_stlo)/2,'degree' );
                azz1(ics) = azimuth( eventcs.evla, eventcs.evlo, sta1_stla,sta1_stlo,'degree' );
                azz2(ics) = azimuth( eventcs.evla, eventcs.evlo, sta2_stla,sta2_stlo,'degree' );
				dazz(ics) = azz1(ics) - azz2(ics);
                rawphv(ics) = ddist(ics)/dtp(ics);
                if (abs(rawphv(ics)-avgphv_evt)/avgphv_evt*100 >= 20 ||  dazz(ics) >= 3.0 )
                    arclen(ics) = -1;
                    rawphv(ics) = -1;
                    isgood(ics) = -1;
                end
                
            else
                arclen(ics) = -1;
                rawphv(ics) = -1;
            end
       
            
       end % end of ics


        %eventraw(ip).phv = rawphv;
        %eventraw(ip).gcarc = arclen;
        %eventraw(ip).isgood = isgood;
        
        goodind = find(isgood>0);
        subplot(3,1,1)
        plot(ddist(goodind),rawphv(goodind),'.','Color',colors(:,ie));
        xlabel('distance between stations pairs (km)');
        ylabel('phv (dtp/ddistance)');
        title(['Period =' ,num2str(periods(ip)),' sec']);
        hold on;
        
        subplot(3,1,2)
        plot(arclen(goodind),rawphv(goodind),'.','Color',colors(:,ie) );
        hold on;
      
        ylabel('phv (dtp/ddistance)');
        xlabel('gcarc (degree)');
        hold on;
        
        subplot(3,1,3)
        plot(arclen(goodind),(rawphv(goodind)-avgphv_evt)/avgphv_evt*100,'.','Color',colors(:,ie));
        ylabel('(phv - evtavgphv)/evtavgphv(%)');
        xlabel('gcarc (degree)');
        hold on;
        
        
    end
    PS = ['T',num2str(periods(ip)),'.ps'];
    print('-dpsc2',PS);
end

