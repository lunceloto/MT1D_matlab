
%反演函数
function [rho,fai]=LMinversion(apprho,eLen,rho,freq,gama,g,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                % 一维大地电磁反演
%                % apprho:观测值（视电阻率）
%                % rho:反演模型电阻率
%                % eLen:反演模型每层层高
%                % freq:计算频率
%                % gama，g:阻尼系数
%                % rho:反演完成后的模型
%                % fai:误差序列
%                % 作者：黄思宁
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nf=length(freq);
nm=length(rho);
A=zeros(nf,nm);

[apprho2, ~]=mt1d_3int(freq,eLen,rho);%计算初始模型的误差
fai(1)=mean(((apprho2-apprho)./apprho).^2);
faimin=fai(1);
nn=0;

if isempty(gcp('nocreate'))
     p=parpool(8);
end

for i=1:n
    [rhoai, ~]=mt1d_3int(freq,eLen,rho);%组装反演矩阵A与B
    B=(apprho'./rhoai'-1);
    drho=0.1;
%     for j=1:nm
%         rho_tmp = rho;                                           %单核版本
%         rho_tmp(j) = rho_tmp(j) * (1 + drho);
%         [rhoci, ~]=mt1d_3int(freq,eLen,rho_tmp);
%         A(:, j) = (rhoci - rhoai) / (rho(j) * drho);
%     end
    parfor j = 1:nm                                           %并行版本
        rho_tmp = rho;
        rho_tmp(j) = rho_tmp(j) * (1 + drho);
        [rhoci, ~] = mt1d_3int(freq, eLen, rho_tmp);
        A(:, j) = (rhoci - rhoai) / (rho(j) * drho);
    end




    X=((A'*A+gama*eye(nm))^(-1))*A'*B;                      %求解迭代量
    X(X>0.5)=0.5;                                         %限制模型修正量，增加收敛稳定性
    X(X<-0.5)=-0.5;
    rho=rho.*(1+X');                                            %修正模型
    [apprho2, ~ ]=mt1d_3int(freq,eLen,rho);         %反演修正模型并求解误差
    fai(i+1)=mean(((apprho2-apprho)./apprho).^2);
    disp([i,fai(i+1)]);


    if fai(i+1)>fai(i)                                      %判断阻尼系数的增减
        gama=gama*g;
    elseif fai(i+1)<fai(i)
        gama=gama/g;
    elseif fai(i+1)==fai(i)
        break;
    end

    if fai(i+1)<0.0001                                       %达到要求的误差时停止
        break;
    end

    figure(3);
    plotmod(rho,eLen,3100);
    title(i);
end

end