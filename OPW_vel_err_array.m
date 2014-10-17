function errs = OPW_vel_err_array(para,localcs)

    v1 = para(1);
    setup_parameters;
    

% in the situation that invert the velocity only
		% gather information
	T = localcs.period;
	stlas = localcs.stlas;
	stlos = localcs.stlos;
	center_sta = localcs.center_sta;
	center_la = stlas(center_sta);
	center_lo = stlos(center_sta);
%     center_la = mean(parameters.lalim);
%     center_lo = mean(parameters.lolim);
	if length(para)==1
		theta1 = azimuth(center_la,center_lo,localcs.evla,localcs.evlo)+180;
	else
		theta1 = para(2);
	end
	% transfer to xy coor
	kx1 = 2*pi/v1/T*cosd(theta1);
	ky1 = 2*pi/v1/T*sind(theta1);
	%stxs = deg2km(stlas-center_la)
	%stys = deg2km((stlos-center_lo)*cosd(center_la))
	[stxs stys] = distance_InArray(center_la,center_lo,stlas, stlos);
    stxs = deg2km(stxs);
    stys = deg2km(stys);
    %calculate the phase misfit
	i = sqrt(-1);
	wave1 = exp(i*(stxs.*kx1+stys.*ky1));
	final_wave = wave1;
	phi_pre = angle(final_wave);
	dphi_pre = phi_pre - phi_pre(center_sta);
	dphi_pre = wrapTo2Pi(dphi_pre);
	dphi_obs = localcs.dtps/T*2*pi;
	dphi_obs = wrapTo2Pi(dphi_obs);
	dph = (dphi_pre - dphi_obs);
	dph = wrapTo2Pi(dph);
	ind = find(dph>pi);
	if ~isempty(ind)
		dph(ind) = dph(ind)-2*pi;
    end
    errs = dph;
    
    badind = find(localcs.isgood<=0);
    errs(badind) = 0;
    
	

%     figure(57)
%     clf
%     hold on
%     plot(localcs.dists,dphi_obs,'x');
%     plot(localcs.dists,dphi_pre,'o');
%     drawnow
%     
%     figure(58)
%     clf
%     hold on
%     plot(stxs,stys,'v');
%     ind = find(abs(errs)>0.04);
%     plot(stxs(ind),stys(ind),'rx');
%     plot(stxs(center_sta),stys(center_sta),'rs');