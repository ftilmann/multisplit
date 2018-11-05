%clc; close all; clear all; clf;
function [data]=do_plot_vars(mcrunno,finalparfi,sacID,filename1,filename2,crted11sac,crted12sac,crted21sac,crted22sac,radialsac,tangensac);
addpath /home/anatolia/UTILS/MatSAC/MatSAC/


!/bin/rm -f tmp.out
!sed 1,1d double_bootstrap.x > tmp.out

data=load('tmp.out');
data1=load(filename1);
data2=load(filename2);

finalpardat=load(finalparfi);

FPDupfinal=finalpardat(:,1);
TDupfinal=finalpardat(:,2);
FPDbotfinal=finalpardat(:,3);
TDbotfinal=finalpardat(:,4);


figure(1)

subplot(4,4,1)
avg1=median(data(:,1)); er1=std(data(:,1)); avgup1=avg1+er1; avglow1=avg1-er1;
hist(data(:,1)); ylabel('Frequency'); xlabel('Top FPDs (s)'); hold on; axis([0 180 0 mcrunno])
vec11=[ones(length([0:0.1:mcrunno]),1).*avg1 [0:0.1:mcrunno]'];
vec12=[ones(length([0:0.1:mcrunno]),1).*avgup1 [0:0.1:mcrunno]'];
vec13=[ones(length([0:0.1:mcrunno]),1).*avglow1 [0:0.1:mcrunno]'];
plot(vec11(:,1),vec11(:,2),'r'); hold on; plot(vec12(:,1),vec12(:,2),'-.r'); plot(vec13(:,1),vec13(:,2),'-.r')


subplot(4,4,5)
plot(data(:,1),data(:,2),'rs'); hold on
plot(data1(:,1),data1(:,2),'+b'); hold on
plot(data2(:,1),data2(:,2),'.k'); hold on
plot(FPDupfinal,TDupfinal,'o','Color',[0 0.2 0.066667]); hold on

ylabel('Top TDs (s)'); xlabel('Top FPDs (deg)');
axis([0 180 0 3])

subplot(4,4,6)
avg2=median(data(:,2)); er2=std(data(:,2)); avgup2=avg2+er2; avglow2=avg2-er2;
hist(data(:,2)); ylabel('Frequency'); xlabel('Top TDs (s)'); hold on; axis([0 3 0 mcrunno])
vec21=[ones(length([0:0.1:mcrunno]),1).*avg2 [0:0.1:mcrunno]'];
vec22=[ones(length([0:0.1:mcrunno]),1).*avgup2 [0:0.1:mcrunno]'];
vec23=[ones(length([0:0.1:mcrunno]),1).*avglow2 [0:0.1:mcrunno]'];
plot(vec21(:,1),vec21(:,2),'r'); hold on; plot(vec22(:,1),vec22(:,2),'-.r'); plot(vec23(:,1),vec23(:,2),'-.r')


subplot(4,4,9)
plot(data(:,1),data(:,3),'rs'); hold on
plot(data1(:,1),data1(:,3),'+b');
plot(data2(:,1),data2(:,3),'.k');
plot(FPDupfinal,FPDbotfinal,'o','Color',[0 0.2 0.066667]); hold on
ylabel('Bottom FPDs (deg)'); xlabel('Top FPDs (deg)');
axis([0 180 0 180])

subplot(4,4,10)
plot(data(:,2),data(:,3),'rs'); hold on
plot(data1(:,2),data1(:,3),'+b');
plot(data2(:,2),data2(:,3),'.k');
plot(TDupfinal,FPDbotfinal,'o','Color',[0 0.2 0.066667]); hold on
ylabel('Bottom FPDs (deg)'); xlabel('Top TDs (s)');
axis([0 3 0 180])

subplot(4,4,11)
avg3=median(data(:,3)); er3=std(data(:,3)); avgup3=avg3+er3; avglow3=avg3-er3;
hist(data(:,3)); ylabel('Frequency'); xlabel('Bottom FPDs (s)'); hold on; axis([0 180 0 mcrunno])
vec31=[ones(length([0:0.1:mcrunno]),1).*avg3 [0:0.1:mcrunno]'];
vec32=[ones(length([0:0.1:mcrunno]),1).*avgup3 [0:0.1:mcrunno]'];
vec33=[ones(length([0:0.1:mcrunno]),1).*avglow3 [0:0.1:mcrunno]'];
plot(vec31(:,1),vec31(:,2),'r'); hold on; plot(vec32(:,1),vec32(:,2),'-.r'); plot(vec33(:,1),vec33(:,2),'-.r')



subplot(4,4,13)
plot(data(:,1),data(:,4),'rs'); hold on
plot(data1(:,1),data1(:,4),'+b');
plot(data2(:,1),data2(:,4),'.k');
plot(FPDupfinal,TDbotfinal,'o','Color',[0 0.2 0.066667]); hold on
ylabel('Bottom TDs (s)'); xlabel('Top FPDs (deg)');
axis([0 180 0 3])

subplot(4,4,14)
plot(data(:,2),data(:,4),'rs'); hold on
plot(data1(:,2),data1(:,4),'+b');
plot(data2(:,2),data2(:,4),'.k');
plot(TDupfinal,TDbotfinal,'o','Color',[0 0.2 0.066667]); hold on
ylabel('Bottom TDs (s)'); xlabel('Top TDs (s)');
axis([0 3 0 3])

subplot(4,4,15)
plot(data(:,3),data(:,4),'rs'); hold on
plot(data1(:,3),data1(:,4),'+b');
plot(data2(:,3),data2(:,4),'.k');
plot(FPDbotfinal,TDbotfinal,'o','Color',[0 0.2 0.066667]); hold on
ylabel('Bottom TDs (s)'); xlabel('Bottom FPDs (deg)');
axis([0 180 0 3])

subplot(4,4,16)
avg4=median(data(:,4)); er4=std(data(:,4)); avgup4=avg4+er4; avglow4=avg4-er4;
hist(data(:,4)); ylabel('Frequency'); xlabel('Bottom TDs (s)'); hold on; axis([0 3 0 mcrunno])
vec41=[ones(length([0:0.1:mcrunno]),1).*avg4 [0:0.1:mcrunno]'];
vec42=[ones(length([0:0.1:mcrunno]),1).*avgup4 [0:0.1:mcrunno]'];
vec43=[ones(length([0:0.1:mcrunno]),1).*avglow4 [0:0.1:mcrunno]'];
plot(vec41(:,1),vec41(:,2),'r'); hold on; plot(vec42(:,1),vec42(:,2),'-.r'); plot(vec43(:,1),vec43(:,2),'-.r')


set(gcf, 'Color','w','Position', [500, 500, 800, 800])


print -dpsc tmp1.ps

figure(2)

subplot(2,2,1)

avg1=median(data(:,1)); er1=std(data(:,1)); avgup1=avg1+er1; avglow1=avg1-er1;
hist(data(:,1)); ylabel('Frequency'); xlabel('Top FPDs (s)'); hold on; axis([0 180 0 1000])
vec11=[ones(length([0:0.1:mcrunno]),1).*avg1 [0:0.1:mcrunno]'];
vec12=[ones(length([0:0.1:mcrunno]),1).*avgup1 [0:0.1:mcrunno]'];
vec13=[ones(length([0:0.1:mcrunno]),1).*avglow1 [0:0.1:mcrunno]'];
plot(vec11(:,1),vec11(:,2),'r'); hold on; plot(vec12(:,1),vec12(:,2),'-.r'); plot(vec13(:,1),vec13(:,2),'-.r')

subplot(2,2,2)

avg2=median(data(:,2)); er2=std(data(:,2)); avgup2=avg2+er2; avglow2=avg2-er2;
hist(data(:,2)); ylabel('Frequency'); xlabel('Top TDs (s)'); hold on; axis([0 3 0 1000])
vec21=[ones(length([0:0.1:mcrunno]),1).*avg2 [0:0.1:mcrunno]'];
vec22=[ones(length([0:0.1:mcrunno]),1).*avgup2 [0:0.1:mcrunno]'];
vec23=[ones(length([0:0.1:mcrunno]),1).*avglow2 [0:0.1:mcrunno]'];
plot(vec21(:,1),vec21(:,2),'r'); hold on; plot(vec22(:,1),vec22(:,2),'-.r'); plot(vec23(:,1),vec23(:,2),'-.r')


subplot(2,2,3)

avg3=median(data(:,3)); er3=std(data(:,3)); avgup3=avg3+er3; avglow3=avg3-er3;
hist(data(:,3)); ylabel('Frequency'); xlabel('Bottom FPDs (s)'); hold on; axis([0 180 0 1000])
vec31=[ones(length([0:0.1:mcrunno]),1).*avg3 [0:0.1:mcrunno]'];
vec32=[ones(length([0:0.1:mcrunno]),1).*avgup3 [0:0.1:mcrunno]'];
vec33=[ones(length([0:0.1:mcrunno]),1).*avglow3 [0:0.1:mcrunno]'];
plot(vec31(:,1),vec31(:,2),'r'); hold on; plot(vec32(:,1),vec32(:,2),'-.r'); plot(vec33(:,1),vec33(:,2),'-.r')

subplot(2,2,4)

avg4=median(data(:,4)); er4=std(data(:,4)); avgup4=avg4+er4; avglow4=avg4-er4;
hist(data(:,4)); ylabel('Frequency'); xlabel('Bottom TDs (s)'); hold on; axis([0 3 0 1000])
vec41=[ones(length([0:0.1:mcrunno]),1).*avg4 [0:0.1:mcrunno]'];
vec42=[ones(length([0:0.1:mcrunno]),1).*avgup4 [0:0.1:mcrunno]'];
vec43=[ones(length([0:0.1:mcrunno]),1).*avglow4 [0:0.1:mcrunno]'];
plot(vec41(:,1),vec41(:,2),'r'); hold on; plot(vec42(:,1),vec42(:,2),'-.r'); plot(vec43(:,1),vec43(:,2),'-.r')

set(gcf,'Color','w')

% plots raw and anisotropy corrected data and their particle motions


[Rtime,R,RSAChdr] = fget_sac(radialsac);
[Ttime,T,TSAChdr] = fget_sac(tangensac);
[C11time,C11,RSAChdr] = fget_sac(crted11sac);
[C12time,C12,TSAChdr] = fget_sac(crted12sac);
[C21time,C21,RSAChdr] = fget_sac(crted21sac);
[C22time,C22,TSAChdr] = fget_sac(crted22sac);




head1=sac(radialsac); t4=head1(3,5); t41=t4 - 15; t42=t4 + 25; 
a1=find(Rtime>=t41 & Rtime<=t42); radtime=Rtime(a1); radial=R(a1);
a2=find(Ttime>=t41 & Ttime<=t42); tantime=Ttime(a2); tangential=T(a2);
a3=find(C11time>=t41 & C11time<=t42); c11time=C11time(a3); c11=C11(a3);
a4=find(C12time>=t41 & C12time<=t42); c12time=C12time(a4); c12=C12(a4);
a5=find(C21time>=t41 & C21time<=t42); c21time=C21time(a5); c21=C21(a5);
a6=find(C22time>=t41 & C22time<=t42); c22time=C22time(a6); c22=C22(a6);

print -dpsc tmp2.ps

figure(3)

subplot(3,2,1)
plot(radtime,radial,'r'); hold on; plot(tantime, tangential,'k'); hold on;
xlabel('time (s)'); legend('radial','tangential'); title(sacID)
subplot(3,2,3)
plot(c11time,c11,'r'); hold on; plot(c12time, c12,'k'); hold on;
xlabel('time (s)'); legend('corrected11','corrected12');
subplot(3,2,5)
plot(c21time,c21,'r'); hold on; plot(c22time, c22,'k'); hold on;
xlabel('time (s)'); legend('corrected21','corrected22'); 

subplot(3,2,2)
plot(radial,tangential)
title('PM')
subplot(3,2,4)
plot(c11,c12)

subplot(3,2,6)
plot(c21,c22)



%sets the size of figure
set(gcf, 'Color','w','Position', [400, 400, 800, 800])


print -dpsc tmp3.ps





