function [ps,phase]=MT1D_Loyar_fwd(rho,h,T)
% The MT 1D forward of the loyar model
% [ps,phase]=MT1D_Loyar_fwd(rho,h,w)
% rho: resistivity
% h: thickness
% fre: Frequency
% ps: Apparenet resistivity
% phase: Impedance phase

t=length(T);
n=length(rho);

Z=zeros(1,n);L=zeros(1,n-1);ps1=zeros(t,1);phase=zeros(t,1);

% initialise the wave-number and the Impedance
k=Z;Z0=Z;

% the dielectric coefficient of air
u0=4*pi*1e-7;
omga=2*pi./T;

for b=1:t
    
    for m=1:n
        k(m)=sqrt(-1i*u0*omga(b)/(rho(m)));
        Z0(m)=-1i*u0*omga(b)/(k(m));
    end
    
    Z(n)=Z0(n);
    while m>1
        L(m)=(Z0(m-1)-Z(m))/(Z0(m-1)+Z(m));
        Z(m-1)=Z0(m-1)*((1-L(m)*exp(-2*k(m-1)*h(m-1)))/(1+L(m)*exp(-2*k(m-1)*h(m-1))));
        m=m-1;
    end
    
    ps1(b)=(abs(Z(1))).^2/(u0*omga(b));
    phase(b)=atan(abs(imag(Z(1))/real(Z(1))))*180/pi;
    
end

ps=ps1;

end
