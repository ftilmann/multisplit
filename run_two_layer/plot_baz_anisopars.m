clear all; close all; clc; clf;

data=load('baz_anisopar.dat');

avgfpdpUP=median(data(:,1));  stdfpdUP=sqrt(sum((((data(:,1)-avgfpdpUP).^2)))/length(data-1)); bound1=avgfpdpUP-stdfpdUP; bound2=avgfpdpUP+stdfpdUP;
avgfpdpLOW=median(data(:,3)); stdfpdLOW=sqrt(sum((((data(:,3)-avgfpdpLOW).^2)))/length(data-1)); bound3=avgfpdpLOW-stdfpdLOW; bound4=avgfpdpLOW+stdfpdLOW;

figure(1)
subplot(2,2,1)
plot([0:1:360],ones(length([0:1:360]),1).*avgfpdpUP,'LineWidth',4,'Color',[.5 .5 .5]); hold on
plot([0:1:360],ones(length([0:1:360]),1).*bound1,'LineWidth',1); hold on
plot([0:1:360],ones(length([0:1:360]),1).*bound2,'LineWidth',1); hold on

scatter(data(:,5),data(:,1),100,data(:,2),'filled','s'); axis([0 360 0 180]); colorbar; 
xlabel('BAZ (deg)'); ylabel('FPD_u_p (deg)')


subplot(2,2,2)

plot([0:1:360],ones(length([0:1:360]),1).*avgfpdpLOW,'LineWidth',4,'Color',[.5 .5 .5]); hold on
plot([0:1:360],ones(length([0:1:360]),1).*bound3,'LineWidth',1); hold on
plot([0:1:360],ones(length([0:1:360]),1).*bound4,'LineWidth',1); hold on
scatter(data(:,5),data(:,3),100,data(:,4),'filled','s'); axis([0 360 0 180]); colorbar;
xlabel('BAZ (deg)'); ylabel('FPD_l_o_w (deg)')


