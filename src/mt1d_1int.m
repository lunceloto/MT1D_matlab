function [apprho,appphs]=mt1d_1int(freq,eLenn,rhon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                % 一维大地电磁一阶有限元正演
%                % apprho:视电阻率
%                % appphs:视相位
%                % freq:计算频率
%                % eLenn:每层层高
%                % rhon:每层电阻率（包含最底层）
%                % 作者：黄思宁
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
miu0= 4 *pi*1e-7;
n=length(rhon)-1;
freq_n=length(freq);

K1=sparse(n+1,n+1);                                                        % 组装求解矩阵K1
for i=1:n
    K1(i:i+1,i:i+1)=K1(i:i+1,i:i+1)+[1,-1;-1,1]/eLenn(i);
end
apprho=zeros(1,freq_n);
appphs=zeros(1,freq_n);
for i=1:freq_n
    omega=2*pi*freq(i);
    K2=sparse(n+1,n+1);
    K3=sparse(n+1,n+1);
    for j=1:n                                   % 组装求解矩阵K2
         K2(j:j+1,j:j+1)=K2(j:j+1,j:j+1)+1i*omega*miu0*eLenn(j)*[1/3,1/6;1/6,1/3]/rhon(j);
    end
    K3 (n+1,n+1)=sqrt(-1i*omega*miu0/rhon(n+1));% 组装求解矩阵K3
    A=K1-K2+K3;
    clear K2 K3;
    A(1,1)=A(1,1)+10^10;
    b=zeros(n+1,1);
    b(1)=A(1,1);
    U=(A^(-1))*b;
    Zimp=1i*omega*miu0*U(1)*eLenn(1)/(U(2)-U(1));%两点微分
    % Zimp=1i*omega*miu0*U(1)/(  -U(1)*( eLenn(2)*eLenn(3)+eLenn(2)*eLenn(4)+eLenn(3)*eLenn(4) ) / ( eLenn(2)*eLenn(3)*eLenn(4) )...
    % +U(2)*( eLenn(3)*eLenn(4) ) / ( eLenn(2)^3-(eLenn(3)+eLenn(4))*eLenn(2)^2+eLenn(2)*eLenn(3)*eLenn(4) )...
    % +U(3)*( eLenn(2)*eLenn(4) ) / ( eLenn(3)^3-(eLenn(2)+eLenn(4))*eLenn(3)^2+eLenn(2)*eLenn(3)*eLenn(4) )...
    % +U(4)*( eLenn(2)*eLenn(3) ) / ( eLenn(4)^3-(eLenn(2)+eLenn(3))*eLenn(4)^2+eLenn(2)*eLenn(3)*eLenn(4) )  );%四点微分
    apprho(i)=(abs(Zimp)^2)/(omega*miu0);
    appphs(i)=atan2(-imag(Zimp), real(Zimp))*180/pi;
    clear A K2 K3 b;
end

end