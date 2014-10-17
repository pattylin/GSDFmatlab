function plot_win_select(event,periods,winpara)
% This function is used to automatically select the window range used for gsdf method.
% The output format is
% v1 = winpara(1); t1 = winpara(2);
% v2 = winpara(3); t2 = winpara(4);
% and the window is defined by L/v1+t1 -- L/v2+t2

if ~exist('periods','var')
    setup_parameters;
    periods = parameters.periods;
end
if ~exist('winpara','var')
    winpara = event.winpara;
end

% Parameters for subplot
Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx;
height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
load('coast.mat');
coast.lat = lat;
coast.lon = long;
clear lat long

freqs = 1./periods;
minf = min(freqs);
maxf = max(freqs);

isdebug = 1;
issaved = 1;

% COLORS!
ORANGE=[255 69 0]./255;

f1=figure(1)
clf
hold on

amp = 0.4;
isgood = [event.stadata(:).isgood];
goodind = find(isgood > 0);
dist = [event.stadata(goodind).dist];
amp = (max(dist)-min(dist))/10*amp;
yrange = [ min(dist)-2*amp max(dist)+2*amp ];
trange = [ mean(dist)/5*0.8  mean(dist)/2*1.2 ];
trange = [ 0  mean(dist)/2*1.2 ];

for ista = 1:length(event.stadata)
    % set up time axis
    if event.stadata(ista).isgood < 0
        continue;
    end
    bgtime = event.stadata(ista).otime - event.otime;
    dt = event.stadata(ista).delta;
    Nt = length(event.stadata(ista).data);
    taxis = bgtime + [0:Nt-1]'*dt;
    data = event.stadata(ista).data;
    fN = 1/2/dt;
    [b,a] = butter(2,[minf/fN, maxf/fN]);
    data = filtfilt(b,a,data);
    data =  data./max(abs(data));
    plot(taxis,data*amp+event.stadata(ista).dist,'-k');
    
    [gausf,faxis] = build_gaus_filter(freqs,dt,Nt,0.06,0.1);
    % get original data and make the fourier transform
    odata = event.stadata(ista).data;
    if size(odata,1) == 1  % in matlab, the fast direction is column
        odata = odata';
    end
    fftodata = fft(odata);
    
    clear envelop_nbands nband nbands norm_envelop
    % apply narrow-band filters
    for ip = 1:length(freqs);
        nband = fftodata .* [gausf(:,ip); zeros(Nt-length(gausf(:,ip)),1)];
        nband = ifft(nband);
        norm_nband = abs(nband)./max(abs(nband));
        %         plot(taxis,norm_nband*amp+event.stadata(ista).dist,'k');
    end % end of loop ip
    stlat = event.stadata(ista).stla;
    stlon = event.stadata(ista).stlo;
end % end of loop sta
if length(winpara)==4
    plot([yrange/winpara(1)+winpara(2)],yrange,'r','linewidth',2);
    disp('R1 xvalues');
    disp([yrange/winpara(1)+winpara(2)])
    disp('R2 xvalues');
    disp([yrange/winpara(3)+winpara(4)])
    disp('R1 and R2 yvalues');
    disp(yrange)
    plot([yrange/winpara(3)+winpara(4)],yrange,'r','linewidth',2);
end

% plot([yrange/5],yrange,'k');
% plot([yrange/2],yrange,'k');
ylim(yrange)
xlim(trange)
% xlim([500 4500]);

title(sprintf('Event: %s Dist: %i ',event.id,mean(dist)));
hold on
drawnow;

%% Draw subplot showing great circle path
[gcplat,gcplon] = track2(event.evla,event.evlo,stlat,stlon);
ix = 1;
iy = 2;
left = sidegap + 1.3*(iy-1)*(vgap+width);
bot = botgap + (ix-1)*(hgap+height);
H = axes('position',[left,bot-.1,width/1.5,height/1.5]);
ax = worldmap([-70 70],[0 360]);
setm(H,'ffacecolor',[1 1 1],'mlabellocation',0,'plabellocation',0);
hold on
plotm(coast.lat,coast.lon,'-b');
plotm(event.evla,event.evlo,'go','markersize',10,'linewidth',2,'markerfacecolor','g');
plotm(stlat,stlon,'rv','markersize',7,'linewidth',2);
plotm(gcplat,gcplon,'-k','linewidth',2);

%% Save the image
if issaved
    print(f1,[num2str(event.id),'.ps'],'-dpsc');
end
