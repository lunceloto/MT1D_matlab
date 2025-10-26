function plotmod(rho,eLen,H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                % eLen:每层层高
%                % rho:每层电阻率（包含最底层）
%                % H:绘制深度（建议比模型总厚度深）
%                % 作者：黄思宁
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn=2000;
eLen=[eLen,H-sum(eLen)];
dz=sum(eLen)/nn;
eLen_n=ceil(eLen/dz);
x=1:dz:sum(eLen);
y=zeros(1,nn);
aaa=0;
for i=1:length(rho)
    y(aaa+1:aaa+eLen_n(i))=rho(i);
    aaa=aaa+eLen_n(i);
end

plot(x(1:nn),y(1:nn));
end