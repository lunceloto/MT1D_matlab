function [Hn,rhon] = make1Dmod(H,rho,min_val, n, total_sum)
% 生成层高按指数递增一维模型，并将原模型电阻率赋值给生成的模型
%       H：          原模型层高
%       rho：        原模型电阻率
%       min_val：    最薄层层高
%       n：          总层数
%       total_sum：  总层高（不含最后一层）
%       Hn：         每层层高
%       rhon：       每层电阻率

    % 参数验证
    if n <= 0 || ~isnumeric(n) || n ~= round(n)
        error('元素个数必须是正整数');
    end
    
    if min_val <= 0
        error('最小值必须大于0');
    end
    
    if total_sum <= min_val * n
        error('总和必须大于最小值乘以元素个数');
    end
    
    % 定义方程：sum(a * r^(k-1)) = total_sum，其中a=min_val
    % 这是一个等比数列求和问题：min_val * (r^n - 1)/(r - 1) = total_sum
    equation = @(r) min_val * (r.^n - 1)./(r - 1) - total_sum;
    
    % 寻找方程的根（r > 1，因为要递增）
    try
        r = fzero(equation, [1 + eps, 100]); % 搜索范围从略大于1到100
    catch
        error('无法找到合适的公比，请检查输入参数');
    end
    
    % 生成指数递增数组
    k = 0:n-1;
    Hn = min_val * r.^k;
    
    % 验证结果（可选，用于调试）
    actual_sum = sum(Hn);
    if abs(actual_sum - total_sum) > 1e-10
        warning('实际总和与目标总和有微小差异: %.10f', abs(actual_sum - total_sum));
    end
    
%     nn=1;
%     Ho=H(1);
%     rhon=zeros(1,n+1);
%     for i=1:n
%         if Ho>=sum(Hn(1:i))
%             rhon(i)=rho(nn);
%         else
%             if nn<length(H)
%                 nn=nn+1;
%                 Ho=sum(H(1:nn));
%             elseif nn<=length(H)+1
%                 Ho=sum(H(1:nn-1));
%             end
%             rhon(i)=rho(nn);
%         end
%     end
%     if Ho>=sum(Hn)
%         rhon(n+1)=rho(nn);
%     else
%         if nn<length(H)
%             nn=nn+1;
%             Ho=sum(H(1:nn));
%         elseif nn<=length(H)+1
%             Ho=sum(H(1:nn-1));
%         end
%         rhon(i)=rho(nn);
%     end
    cumH  = cumsum(H);
    cumHn = cumsum(Hn);

    nn = 1;
    rhon = zeros(1, n + 1);

    for i = 1:n
        if cumH(nn) >= cumHn(i)
            rhon(i) = rho(nn);
        else
            rhon(i) = rho(nn+1);
            nn = min(nn + 1, length(H)); % 防止越界
        end
    end

    % 处理最后一个点
    if cumH(nn) >= cumHn(end)
        rhon(n+1) = rho(nn);
    else
        rhon(n+1) = rho(nn+1);
    end
end