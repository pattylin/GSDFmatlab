%% This script is similar to stack_phv.m, but adding azimuthal anisotropy inversion into the result.
%
% written by Ge Jin, LDEO
% jinwar@gmail.com, ge.jin@ldeo.columbia.edu
%


clear;

phase_v_path = './eikonal/'
r = 0.03;
isfigure = 1;

setup_parameters

comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
min_csgoodratio = parameters.min_csgoodratio;
min_phv_tol = parameters.min_phv_tol;
max_phv_tol = parameters.max_phv_tol;
is_raydense_weight = parameters.is_raydense_weight;
min_event_num = parameters.min_event_num;
err_std_tol = parameters.err_std_tol;
smsize = parameters.smsize;
off_azi_tol = parameters.off_azi_tol;
is_one_phi = parameters.is_one_phi;

xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);

[xi yi]=ndgrid(xnode,ynode);

for ip=1:length(periods)
	avgphv(ip).sumV = zeros(Nx,Ny);
	avgphv(ip).sumweight = zeros(Nx,Ny);
	avgphv(ip).GV_std = zeros(Nx,Ny);
	avgphv(ip).eventnum = zeros(Nx,Ny);
	avgphv(ip).xi = xi;
	avgphv(ip).yi = yi;
	avgphv(ip).xnode = xnode;
	avgphv(ip).ynode = ynode;
	avgphv(ip).period = periods(ip);
end


% read in bad event list, if existed
if exist('badevt.lst')
    badevts = textread('badevt.lst','%s');
    disp('Found Bad event')
    disp(badevts)
end


phvmatfiles = dir([phase_v_path,'/*_',comp,'.mat']);
eventnum =  length(phvmatfiles);

GV_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
azi_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));
raydense_mat = zeros(Nx,Ny,length(phvmatfiles),length(periods));

% Gather information
for ie = 1:length(phvmatfiles)
	if ~exist('badevt.lst') || isempty(find(strncmpi(phvmatfiles(ie).name,badevts,12), 1))
        temp = load([phase_v_path,phvmatfiles(ie).name]);
        eventphv = temp.eventphv;
        disp(['ie: ',num2str(ie), ' event:' ,eventphv(1).id]);
        evla(ie) = eventphv(ip).evla;
        evlo(ie) = eventphv(ip).evlo;
        for ip=1:length(periods)
            GV_mat(:,:,ie,ip) = eventphv(ip).GV;
            raydense_mat(:,:,ie,ip) = eventphv(ip).raydense;
            azi = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
            azi = rad2deg(azi);
            azi_mat(:,:,ie,ip) = azi;
        end
        gcazi_mat(:,:,ie) = azimuth(xi,yi,evla(ie),evlo(ie));
    else
        phvmatfiles(ie).name
        GV_mat(:,:,ie,:) = NaN;
        raydense_mat(:,:,ie,:) = NaN;
        azi_mat(:,:,ie,:) = NaN;
        gcazi_mat(:,:,ie) = NaN;
    end
end
figure(11)
clf

for ip = 1:length(periods)
    isophv=zeros(Nx,Ny);
    isophv_std=zeros(Nx,Ny);
    aniso_strength=zeros(Nx,Ny);
    aniso_strength_std=zeros(Nx,Ny);
    aniso_azi=zeros(Nx,Ny);
    aniso_azi_std=zeros(Nx,Ny);
    aniso_1phi_strength=zeros(Nx,Ny);
    aniso_1phi_azi=zeros(Nx,Ny);	
	% start to inverse azimuthal anisotropy grid by grid
	for mi=1:Nx
		disp(['ip:',num2str(ip),' process: ',num2str(mi/Nx)]);
        for mj=1:Ny
            n=0;
            clear phV_best azi phV dist gcazi
			for ie = 1:eventnum
				avgV=GV_mat(mi,mj,ie,ip);
                lowi=max(1,mi-smsize);
                upi=min(Nx,mi+smsize);
                lowj=max(1,mj-smsize);
                upj=min(Ny,mj+smsize);
                for ii=lowi:upi
                    for jj=lowj:upj
                        if ~isnan(GV_mat(ii,jj,ie,ip))
                            n=n+1;
                            azi(n)=azi_mat(ii,jj,ie,ip);
							gcazi(n) = gcazi_mat(ii,jj,ie);
                            phV(n)=GV_mat(ii,jj,ie,ip);                            
                            %debug
							%if n==15 || n== 20
							%	disp(['ievent=',num2str(ie)]);
							%end
                        end
                    end
                end

                
			end  % enent loop

			% needed  pylin.patty 2014.04 ------------------------ 
            if n <= min_event_num*((2*smsize).^2)
                isophv(mi,mj)=NaN;
                isophv_std(mi,mj)=NaN;
                aniso_strength(mi,mj)=NaN;
                aniso_strength_std(mi,mj)=NaN;
                aniso_azi(mi,mj)=NaN;
                aniso_azi_std(mi,mj)=NaN;
                continue;
            end 
			% ----------------------------------------------------
			
			% Get rid of significant off circle events
			diffazi = azi - gcazi;
			ind = find(diffazi>180);
			if ~isempty(ind)
				diffazi(ind) = diffazi(ind) - 360;
			end
			ind = find(diffazi<-180);
			if ~isempty(ind)
				diffazi(ind) = diffazi(ind) + 360;
			end
			ind = find(abs(diffazi)> off_azi_tol);
            
			if ~isempty(ind)
				azi(ind) = [];
				phV(ind) = [];
            end
            % add pylin.patty 2014.04 -------
			% if phV are far away from the mean by 5%, remove those points before fitting azi anisotropy
            ind = find(abs(phV-mean(phV))/mean(phV)*100 >=5 );
			disp(['ip :',num2str(ip),'    mi,mj: ',num2str(mi),',',num2str(mj), '  ind: ',num2str(ind)]);
            %ind = find(abs(phV-mean(phV)) >= 0.05*mean(phV) )
             
			%ind = find(abs(phV-mean(phV)) > err_std_tol*std(phV) )
  
            if ~isempty(ind)
				azi(ind) = [];
				phV(ind) = [];
            end
            % ------------------------------
            
 
            if is_one_phi && length(azi) < 5 
				isophv(mi,mj)=NaN;
                isophv_std(mi,mj)=NaN;
                aniso_strength(mi,mj)=NaN;
                aniso_strength_std(mi,mj)=NaN;
                aniso_azi(mi,mj)=NaN;
                aniso_azi_std(mi,mj)=NaN;
                numer_phv(mi,mj) = NaN; 
                continue;
			else
            	if is_one_phi
                	[para fiterr]=fit_azi_anisotropy_1phi(azi,phV);
            	else
               	 [para fiterr]=fit_azi_anisotropy(azi,phV);
				end
				numer_phv(mi,mj) = length(azi);
				parastd=confint(para);
				isophv(mi,mj)=para.a;
				isophv_std(mi,mj)=parastd(2,1)-parastd(1,1);
				aniso_strength(mi,mj)=para.d;
				aniso_azi(mi,mj)=para.e;
				if para.e > 180
					aniso_azi(mi,mj)=para.e-180;
				elseif para.e < 0
					aniso_azi(mi,mj)=para.e+180;
				end
				if is_one_phi
					aniso_1phi_strength(mi,mj)=para.b;
					aniso_1phi_azi(mi,mj)=para.c;
					aniso_strength_std(mi,mj)=parastd(2,4)-parastd(1,4);
					aniso_azi_std(mi,mj)=parastd(2,5)-parastd(1,5);
				else
					aniso_strength_std(mi,mj)=parastd(2,2)-parastd(1,2);
					aniso_azi_std(mi,mj)=parastd(2,3)-parastd(1,3);
				end          
			
				if isfigure && mi == 2 && mj == 2 
					plot_position = [ 1 3 5 7 2 4 6 8 ];
					subplot(4,2,plot_position(ip))
					hold on
					allazi = -200:200;
					if ~is_one_phi 
						plot(allazi,(para.a*(1+para.d*cosd(2*(allazi-para.e)))-mean(phV))/mean(phV)*100,'r','linewidth',2)
					elseif is_one_phi 
						plot(allazi,(para.a*(1+para.b*cosd(allazi-para.c)+para.d*cosd(2*(allazi-para.e)))-mean(phV))/mean(phV)*100,'r','linewidth',2)
					end
					plot(azi,(phV-mean(phV))/mean(phV)*100,'x','linewidth',2);
					%length(azi)
					max(abs((phV-mean(phV))/mean(phV)*100))
					%[azi (phV-mean(phV))/mean(phV)*100]
					xlim([-180 180])
					set(gca,'xtick',[-180 -135 -90 -45 0 45 90 135 180])
					title(['Period= ',num2str(periods(ip)),' (s)'])
					ylim([-6 6])
					set(gca,'ytick',[-4:2:4])
					%pause;
					hold on;
			    end
			end
		end  % mj loop
	end % mi loop
	avgphv_aniso(ip).isophv = isophv;
	avgphv_aniso(ip).isophv_std = isophv_std;
	avgphv_aniso(ip).aniso_strength = aniso_strength;
	avgphv_aniso(ip).aniso_azi = aniso_azi;
	avgphv_aniso(ip).parameters = parameters;
	avgphv_aniso(ip).xi = xi;
	avgphv_aniso(ip).yi = yi;
	avgphv_aniso(ip).numer_phvmeasurment = numer_phv;

	if is_one_phi
		avgphv_aniso(ip).aniso_1phi_strength=aniso_1phi_strength;
		avgphv_aniso(ip).aniso_1phi_azi=aniso_1phi_azi;
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	else
		avgphv_aniso(ip).aniso_strength_std=aniso_strength_std;
		avgphv_aniso(ip).aniso_azi_std=aniso_azi_std;
	end          
end % end of period loop
                	psfile =  ['azimuthal_cosine.ps']
                	print('-dpsc2',psfile)

% save to output filename
filename = ['eikonal_stack_aniso_',comp,'.mat'];
save(filename,'avgphv_aniso');


N=3; M = floor(length(periods)/N)+1;

% === plot avgphv_aniso.isophv  =================================================
figure(56)
clf
for ip = 1:length(periods)
	subplot(M,N,ip);
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
	h1=surfacem(xi,yi,avgphv_aniso(ip).isophv);
	colorbar
	load seiscmap
	colormap(seiscmap)
	avgv = nanmean(avgphv_aniso(ip).isophv(:));
	caxis([avgv*(1-r) avgv*(1+r)])
end

% === plot avgphv_aniso.strength  =================================================
figure(57)
clf
for ip = 1:length(periods)
	subplot(M,N,ip);
	ax = worldmap(lalim, lolim);
	set(ax, 'Visible', 'off')
%	h1=surfacem(xi,yi,avgphv_aniso(ip).isophv);
 	h1=surfacem(xi,yi,avgphv_aniso(ip).aniso_strength);
	colorbar
	load seiscmap
	colormap(seiscmap)
% 	avgv = nanmean(avgphv_aniso(ip).isophv(:));
%	caxis([avgv*(1-r) avgv*(1+r)])
     caxis([0 0.05]);
	u=avgphv_aniso(ip).aniso_strength.*cosd(avgphv_aniso(ip).aniso_azi)*50;
	v=avgphv_aniso(ip).aniso_strength.*sind(avgphv_aniso(ip).aniso_azi)*50./cosd(mean(lalim));
	[m n]=size(xi);
	for ix=1:m
		for iy=1:n
			if avgphv_aniso(ip).aniso_azi_std(ix,iy) < 40 && avgphv_aniso(ip).aniso_strength(ix,iy)>0.01
				h=plotm([xi(ix,iy)-u(ix,iy)/2 xi(ix,iy)+u(ix,iy)/2],...
					[yi(ix,iy)-v(ix,iy)/2 yi(ix,iy)+v(ix,iy)/2],'k-');
 				set(h,'linewidth',2)
			end
		end
	end
end

