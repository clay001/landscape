%cdf x is the time series, time*genes; i is the ith gene; 
% 函数mycdfDiff
% 输入（7000*48的模拟结果,438*48的真实表达矩阵,48*7的真实轨迹,fit的高斯模型参数）
% 输出（1*48的cdfDiff）
function diffCDF=mycdfDiff(xSimulate,xReal,realTraj,fittingDataTemp)
    % 值传递，真实表达矩阵438*48
    xRealNorm=xReal;
    % stage的个数 7 
    numPoints=size(realTraj,2);
    

    %%（1，48）
    diffCDF=zeros(1,size(xRealNorm,2));

    % 对每个矩阵
    for genei=1:size(xRealNorm,2)
        % 这里命名不一样了有点容易搞错
        % 1*48
        % 我觉得这里的range，也就是bin的设计有点问题
        temp=xSimulate(genei,:)';
            
        % hist
        % Y是排序后的temp,I1是元素再原tmp中的index
        [Y,I1]=sort(temp);
        % 再排序从小到大，index是原I1中的index
        [Y2,I2]=sort(I1);
        
        %[N2,X]=hist(xReal(:,genei),Y);
        %DF_real=cumsum(N2/sum(N2));

        % 1*4 第一行是参数
        p = fittingDataTemp((genei-1)*3+1,:);
        % 有效长度
        len= 4-length(find(p==0));
        
        % 第二行是平均值，是0的参数列就不看了，假设有效长度是2，则mu为2*1
        MU = fittingDataTemp((genei-1)*3+2,1:len)';

        % 第三列是sigma值
        if len==1
            % cat沿指定维度串联数组，3维度就是平行增设的意思，按照每一列平行增设
            SIGMA = cat(3,fittingDataTemp((genei-1)*3+3,1));
        elseif len==2
            % 应该是1*1*2
            SIGMA = cat(3,fittingDataTemp((genei-1)*3+3,1),fittingDataTemp((genei-1)*3+3,2));
        elseif len==3
            SIGMA = cat(3,fittingDataTemp((genei-1)*3+3,1),fittingDataTemp((genei-1)*3+3,2),fittingDataTemp((genei-1)*3+3,3));
        elseif len==4
            SIGMA = cat(3,fittingDataTemp((genei-1)*3+3,1),fittingDataTemp((genei-1)*3+3,2),fittingDataTemp((genei-1)*3+3,3),fittingDataTemp((genei-1)*3+3,4));
        end

        % p的参数也去掉0列
        p = fittingDataTemp((genei-1)*3+1,1:len)';
        % 还原生成高斯混合模型
        obj = gmdistribution(MU,SIGMA,p);
        
        % Cumulative distribution function 累计分布函数
        % 喂入高斯混合模型的分布，还有对应评估的点，也就是Y（排序好的数据）
        % 返回一个和Y同形的累计分布函数  48*1
        DF_real= cdf(obj,Y);
        
        % 不推荐，可以改用histogram
        % 直方图，喂入7000个细胞在第一个基因处的向量，Y作为中心bin的坐标
        % 返回的是一个1*48的double，记录了每个值为中心的bin中有多少细胞
        N1=hist(xSimulate(:,genei),Y);
        
        % 归一化，并且做个累加
        DF=cumsum(N1/sum(N1));
        % 所有基因的分布情况和真实分布计算误差平方矩阵
        diffs=(DF_real'-DF).^2;

        % 对应位置填入该基因的误差平方和
        diffCDF(genei)=sum(diffs);
    end
    
    % 是一个1*48的double，描述了cdf的误差
    % 1*48 最后再除以基因的个数和细胞阶段数 （48*7）
    % 归一化也有疑问
    diffCDF=diffCDF/size(xRealNorm,2)/numPoints;
    

end