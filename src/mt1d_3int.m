function [apprho,appphs]=mt1d_3int(freq,eLenn,rhon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                % 一维大地电磁三阶有限元正演
%                % apprho:视电阻率
%                % appphs:视相位
%                % freq:计算频率
%                % eLenn:每层层高
%                % rhon:每层电阻率（包含最底层）
%                % 作者：黄思宁
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mui0= 4 *pi*1e-7;
n=length(rhon)-1;
freq_n=length(freq);

K1=sparse(3*n+1,3*n+1);
for i=1:n                                  % 组装求解矩阵K1
    m=3*i-2;
    K1(m:m+3,m:m+3)=K1(m:m+3,m:m+3)+...
        [148 -189 54 -13;
        -189 432 -297 54;
        54 -297 432 -189;
        -13 54 -189 148]/( 40*eLenn(i) );
    rho_n(m:m+2,1)=rhon(i);
end
rho_n(3*n+1,1)=rho_n(3*n,1);

%%% 求解矩阵K2的系数矩阵
kMat = [357  216  -81    20;    %K11
        216  324  -162   18;    %K12
        -81  -162 81     18;    %K13
        20   18   18     20;    %K14
        216  324  -162   18;    %K21
        324  2187 0      81;    %K22
        -162 0    0      -162;  %K23
        18   81   -162   -81;   %K24
        -81  -162 81     18;    %K31
        -162 0    0      -162;  %K32
        81   0    2187   324;   %K33
        18   -162 324    216;   %K34
        20   18   18     20;    %K41
        18   81   -162   -81;   %K42
        18   -162 324    216;   %K43
        20   -81  216    357;]; %K44
apprho=zeros(1,freq_n);
appphs=zeros(1,freq_n);
for i=1:freq_n
    omega=2*pi*freq(i);
    K2=sparse(3*n+1,3*n+1);
    K3=sparse(n+1,n+1);
    for j=1:n                                 % 组装求解矩阵K2
        m=3*j-2;
        KK=1i*omega*mui0*eLenn(j)*kMat*(1./rho_n(m:m+3))/6720;
        K2(m:m+3,m:m+3)=K2(m:m+3,m:m+3)+reshape(KK,4,4);
    end
    K3(3*n+1,3*n+1)=sqrt(-1i*omega*mui0/rhon(n+1));
    A=K1-K2+K3;
    clear K2 K3;
    A(1,1)=A(1,1)+10^10;
    b=zeros(3*n+1,1);
    b(1)=A(1,1);
    U=(A^(-1))*b;
    Zimp=1i*omega*mui0*U(1)*2*eLenn(1)/(-11*U(1)+18*U(2)-9*U(3)+2*U(4));   %四点微分
    apprho(i)=(abs(Zimp)^2)/(omega*mui0);
    %appphs(i)=atan2 (-imag(Zimp), real(Zimp))*180/pi;
    appphs(i)=atan(abs(imag(Zimp)/real(Zimp)))*180/pi;
    clear A K2 K3 b; 
end

end