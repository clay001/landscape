% x is a vecter
% F函数
% 输入（单个细胞基因表达情况48*1，整体基因的均值和标准差2*48）
% 输出（激活函数作用后的加权值1*48）
function fx=F(x,fittingData)
    if size(x,1)==1
        x=x';
    end
    
    
    %old
    fx = 1./(1+exp(-(x'-fittingData(1,:))./fittingData(2,:)));
    %new    
%     w=fittingDataTemp(1:3:end,:);
%     mu=fittingDataTemp(2:3:end,:);
%     sigma=fittingDataTemp(3:3:end,:);
%     
%     fx = sum(w./(1+exp(-(repmat(x,1,4)-mu)./sigma)),2);
end