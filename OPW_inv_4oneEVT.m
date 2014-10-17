function [sta_phv sta_theta sta_gcarc] = OPW_inv_4oneEVT(EVT,ip)
%function [sta_phv sta_theta sta_gcarc] = OPW_inv_4oneEVT(EVTCSm,ip)

setup_parameters

comp = parameters.component;
load(['./CSmeasure/',num2str(EVT),'_cs_',comp,'.mat']);


load(['./eikonal/',num2str(EVT),'_eikonal_',comp,'.mat']);

%load(EVTCSm);

sample_range = parameters.maxstadist;

evla = eventcs.evla;
evlo = eventcs.evlo;
stlas = eventcs.stlas;
stlos = eventcs.stlos;

opts = optimset('Display','off');

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:');
	disp(badstnms);
end

if exist('badstnms','var')
    badstaids = find(ismember(eventcs.stnms,badstnms));
    goodstaids = find(~ismember(eventcs.stnms,badstnms));
else
    badstaids = [];
    goodstaids = 1:length(stlas);
end


for ista = 9 : 9 %length(stlas)

    if ~ismember(ista,badstaids) 
        eventcs.stnms(ista)
        center_sta = ista;
        center_stla = stlas(ista);
        center_stlo = stlos(ista);

        arclen = distance( evla, evlo, center_stla, center_stlo ,'degree' );


        stadists = distance(stlas,stlos,center_stla,center_stlo);
        stadists = deg2km(stadists);

        nbstaids = find(stadists < sample_range);

        % reconstruct the data format to local array index
        CSnum = 0;
        localcs = [];
        localcs.stlas = stlas(nbstaids);
        localcs.stlos = stlos(nbstaids);
        localcs.dists = deg2km(distance(localcs.stlas,localcs.stlos,evla,evlo));
        localcs.evla = evla;
        localcs.evlo = evlo;
        localcs.center_sta = find(nbstaids == center_sta);
        localcs.period = parameters.periods(ip);

        localcs.isgood(localcs.center_sta) = 1;    

        for jsta = 1:length(nbstaids)
            localcs.dtps(jsta) = 0;
            if ( ~isnan(eventphv(ip).traveltime(jsta)))
                
                localcs.isgood(jsta) = 1;
                localcs.dtps(jsta) = eventphv(ip).traveltime(jsta)-eventphv(ip).traveltime(ista);
            else
                localcs.isgood(jsta) = 0;
                localcs.dtps(jsta) = [];
            end
        end


        
        

        
%         % find the CS measurements that contains the center station
%         for ics = 1:length(eventcs.CS)
%             if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
%                 %ics
%             else
%     
%                 if eventcs.CS(ics).sta1==center_sta 
%                     staid = find(nbstaids==eventcs.CS(ics).sta2);
%                     localcs.dtps(staid) = -eventcs.CS(ics).dtp(ip);
%                     localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
%                 elseif eventcs.CS(ics).sta2==center_sta 
%                     staid = find(nbstaids==eventcs.CS(ics).sta1);
%                     localcs.dtps(staid) = eventcs.CS(ics).dtp(ip);
%                     localcs.isgood(staid) = eventcs.CS(ics).isgood(ip);
%                 end
%             end
%         end % ics

        % define the initial model
        center_la = localcs.stlas(localcs.center_sta);
        center_lo = localcs.stlos(localcs.center_sta);
        [epi_dist baz] = distance(center_la,center_lo,evla,evlo);


        if length(find(localcs.isgood)) > 5
            para0(1) = parameters.refphv(ip);

            para0(2) = baz+180;
            paraL = []; %[para0(1)*0.9 baz+180-0.01];
            paraU = []; %[para0(1)*1.1 baz+180+0.01];
            [para,resnorm,residual,exitflag] = lsqnonlin(@(para) OPW_vel_err_array(para,localcs),para0,paraL,paraU,opts);
            sta_phv(ista) = para(1);
            sta_theta(ista) = para(2);
             
%             temp = polyfit(localcs.dists,localcs.dtps,1);
%             sta_phv(ista) = 1./ temp(1);
%             sta_theta(ista) = para0(2);        
%             figure(56)
%             hold on
%             plot(length(localcs.isgood),para(1),'x');

            if ( sta_theta(ista) > 360 ) 
                sta_theta(ista) = sta_theta(ista) - 360;
            end
            sta_gcarc(ista) = arclen;
        else
            sta_phv(ista) = NaN;
            sta_theta(ista) = NaN;
            sta_gcarc(ista) = NaN;
        end
    else
        sta_phv(ista) = NaN;
        sta_theta(ista) = NaN;
        sta_gcarc(ista) = NaN;   
    end


end
