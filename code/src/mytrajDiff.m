%cdf x is the time series, time*genes; i is the ith gene; 
% 函数mytrajDiff
% 输入（模拟的表达情况7*48，真实表达情况48*7，这个weight是每个基因在不同细胞阶段表达水平的相近程度）
% 输出（每个基因的相对误差和）
function diffTraj=mytrajDiff(xSimulate,realTraj,weight)
    % 细胞状态的数目
    numPoints=size(realTraj,2);
    % 度量模拟数据和真实数据的轨迹差别
    % 前一项是除以细胞状态数目，基因数目（7*48）
    % 后一项是差异加权以后按行求和  48*1
    diffTraj=1/numPoints*1/size(realTraj,1)*sum((xSimulate'-realTraj).^2.*weight,2);

end