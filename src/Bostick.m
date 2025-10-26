function [rho,eLen]=Bostick(rhoT,phsT,freq)
% 大地电磁Bostick反演代码
% 个人学习用
miu=4 *pi*1e-7;
omiga=2*pi*freq;
for i=1:length(rhoT)
    eLen(i)=sqrt(abs(rhoT(i))/(miu*omiga(i)));
    rho(i)=abs(rhoT(i)*(180/(2*phsT(i))-1));
end