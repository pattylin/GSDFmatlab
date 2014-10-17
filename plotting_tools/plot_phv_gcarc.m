% run plot_phv_DC first 
% plot one period, all event evtavgphv v.s. gcarc / baz 


figure(1)
clf
iperiod = 4 
for ie = 1:numbers_events
gcarc(ie) = distance(eventinfo(ie).evla, eventinfo(ie).evlo, 9, -146);
gcabaz(ie) = azimuth(9,-146,eventinfo(ie).evla,eventinfo(ie).evlo);
evdp(ie) = eventinfo(ie).evdp;
end 

subplot(3,1,1)
plot(gcarc, evtavgphv(iperiod).phv,'xb','linewidth',2);
hold on;
title({['Periods: ',num2str(periods(iperiod))]; 'phv vs gcarc'})


periods(iperiod)
%load pa5phvR
%[temp ipt] = min(abs(2*pi./w - periods(iperiod)));
%plot(get(gca,'xlim'), [phv(ipt) phv(ipt)], '--r','linewidth',1);
%plot(get(gca,'xlim'), [sumphv(iperiod).phv  sumphv(iperiod).phv], '--b','linewidth',1);



subplot(3,1,2)
plot(gcabaz, evtavgphv(iperiod).phv,'xb','linewidth',2);
hold on;

title('phv vs baz')
periods(iperiod)
%load pa5phvR
%[temp ipt] = min(abs(2*pi./w - periods(iperiod)));
%plot(get(gca,'xlim'), [phv(ipt) phv(ipt)], '--r','linewidth',1);
%plot(get(gca,'xlim'), [sumphv(iperiod).phv  sumphv(iperiod).phv], '--b','linewidth',1);




subplot(3,1,3)
plot(evdp, evtavgphv(iperiod).phv,'xb','linewidth',2);
hold on;

title('phv vs evdp')
periods(iperiod)
%load pa5phvR
%[temp ipt] = min(abs(2*pi./w - periods(iperiod)));
%plot(get(gca,'xlim'), [phv(ipt) phv(ipt)], '--r','linewidth',1);
%plot(get(gca,'xlim'), [sumphv(iperiod).phv  sumphv(iperiod).phv], '--b','linewidth',1);

