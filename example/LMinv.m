clc; clear; close all;
addpath('..\src');
nf=16;
nm=50;
nN=40;
H=3000;
alpha=0.6;
gama=8;
g=2;
%%%初始模型
eLen0=[500,100,500];
rho0=[100,10,2000,500];
ccc=-4:8/nf:4;
freq=flip(10.^ccc);
%%%正演
[apprho,appphs]=MT1D_Loyar_fwd(rho0,eLen0,flip(freq));

%%%添加10%高斯噪声
apprho=awgn(apprho',10*log10(4));
%save rhoobs apprho appphs;
%load rhoobs.mat;

%%%初始模型
[eLenn,rhon] = make1Dmod(10,[100,100],2, 30, 1500);
%eLenn=(4000/nm)*ones(1,nm);
%rhon=200*ones(1,nm+1);

%%%Bostic反演
%[rhon,eLen]=Bostick(apprho,appphs,freq);
%figure(1);
%plotmod(rho0,eLen0,3100);
%hold on;
%plotmod(rhon,eLen,3100);

%%迭代反演
tic
[rhon,fai]=LMinversion(apprho,eLenn,rhon,freq,gama,g,nN);
toc

figure(1);
plotmod(rho0,eLen0,5100);
hold on;
plotmod(rhon,eLenn,5100);
xlabel('深度');
ylabel('电阻率');
title('反演结果对比图');
legend('正演模型','反演结果')


figure(2)
loglog(fai);
xlabel('迭代次数');
ylabel('误差');
title('误差收敛曲线');


figure(3);
[apprho2,appphs2]=mt1d_3int(freq,eLenn,rhon);
subplot(2,1,1);
loglog(1./freq,apprho);
hold on;
loglog(1./freq,apprho2);
subplot(2,1,2);
loglog(1./freq,appphs);
hold on;
loglog(1./freq,appphs2);



