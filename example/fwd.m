clc;clear;close all;
addpath('..\src');
eLen=[500,100,500];
rho=[100,10,1000,500];
n=40;
H=3000;
alpha=0.6;
ccc=-3:6/n:3;
freq=10.^ccc;
[zn,rhon] = make1Dmod(eLen,rho,1, 50, 1500);
figure(1)
plotmod(rho,eLen,H)
[apprho,appphs]=MT1D_Loyar_fwd(rho,eLen,1./freq);
[apprho1,appphs1]=mt1d_1int(freq,zn,rhon);
[apprho3,appphs3]=mt1d_3int(freq,zn,rhon);

figure(2)
loglog(1./freq,apprho);
hold on;
loglog(1./freq,apprho1,'*');
loglog(1./freq,apprho3,'*');
legend('解析解','一次插值','三次插值');
set(gca,'XTick',[0.00100000000000000	0.0100000000000000	0.100000000000000	1	10	100	1000]);
xlabel('周期/s');
ylabel('视电阻率 ρ');
title('视电阻率对比图');