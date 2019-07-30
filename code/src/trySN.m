%函数trySN
%输入（初始参数，随机生成的起始数据，每个基因的均值和标准差，包含数据的结构体，基因轨迹和权重）
%输出（总误差，7000*48的模拟矩阵，1*1000的模拟cell，7*48的平均模拟矩阵）
function [g,xSimulate1,xSimulate,xSimulate2] = trySN(paramInit,randomXInits,fittingData,hopland,realTraj,weight)
    %%
    orgiFitData=hopland.orgiFitData;
    developLine=hopland.developLine;

    %参数传递
    xReal=orgiFitData;
    %创建全0权重矩阵  
    weightMatrix=zeros(size(xReal,2),size(xReal,2));
    count=1;

    %逐列填充，填充上半部分，用count从初始化的参数列表中取值
    for j=1:size(xReal,2)
        for i=1:j
            weightMatrix(i,j)=paramInit(count);
            % 下半部分设成对称
            weightMatrix(j,i)=weightMatrix(i,j);
            % 往后计数
            count=count+1;
        end
    end
    
    % 看起来很长，其实就是往后再取一倍的基因个数的参数作为sigmaW（48*1）
    sigmaW = paramInit(((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+1):((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+size(xReal,2))); %
    % 
    Ttime=1;
    % 取a
    a=paramInit( ((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+size(xReal,2)+1):((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+2*size(xReal,2)));
    % 取sigma
    sigma=paramInit( ((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+size(xReal,2)*2+1):((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+3*size(xReal,2)));
    
    %% 如果a的长度为1，意味着只有一个基因，不懂为什么分开判断
    if length(a)==1
        A = a*ones(1,size(xReal,2));        
    else
        A = a;       
    end

    % 依然是一个为一判断
    if length(sigmaW)==1
        tempsigmaW = sigmaW*ones(size(xReal,2),size(xReal,2));
    % 若sigmaW的第一个维度为1
    elseif size(sigmaW,1)==1
        % 变为竖排列然后横向复制  genesize * genesize
        tempsigmaW = repmat(sigmaW',1,size(xReal,2));    
    % 否则就原始状态横向复制
    else
        tempsigmaW = repmat(sigmaW,1,size(xReal,2));     
    end
     
    % 同理
    if length(sigma)==1
        sigma = sigma*ones(1,size(xReal,2));
    else
        sigma = sigma;       
    end
    
    % 大小为1*初始细胞个数
    xSimulate=cell(1,size(randomXInits,1));

    % 循环
    for timepoints=1:size(randomXInits,1)
        % 逐个取出初始细胞
        Xzero = randomXInits(timepoints,:);
        % 这里Ttime根据上文是设为1，放进HCN中进行模拟，结果放入Xsimulate中，当作这一个细胞演化成其他cell stage的情况
        xSimulate{timepoints}=hopfieldNetworkContinuousModel(weightMatrix, Ttime, sigma, A, tempsigmaW, Xzero, fittingData, developLine);
    end

    % 转置一下并且转换为mat格式，表征这一堆起始细胞各自的完整成长情况，按初始细胞顺序竖着排列
    % 7000*48
    xSimulate1=cell2mat(xSimulate');
    
    % 大小是stage个数*基因数 7*48
    xSimulate2=zeros(size(xSimulate{1}));
    % 还是按初始细胞个数依次循环
    for timepoints=1:size(randomXInits,1)
        % 对这1000次独立的模拟加和，形状仍然是7*48
        xSimulate2=xSimulate2+xSimulate{timepoints};
    end

    % 除以1000，求平均，7*48
    xSimulate2=xSimulate2/size(randomXInits,1);
    

    %%得到每个基因的加和差距
    x1=mytrajDiff(xSimulate2,realTraj,weight);
    
%     x2=mycdfDiff(xSimulate1,xReal,realTraj,fittingDataTemp);

    % 再得到总误差
    g=sum(x1);
 

end

