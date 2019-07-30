%--------------------------------------------------
% HopLand algorithm
%--------------------------------------------------

function hopland=runHopLand(hopland,givenStartPoint,ifdoComparison)
    %% Step 1: fit muliGuassian
    %得到每个基因的均值标准差，fit的高斯模型的参数
    [fittingData,fittingDataTemp]=fitMixtureGaussian(hopland);
    
    %% TRAINING
    %开始训练，如果时间信息比较丰富
    if hopland.ifTimeseries        
        % Step 2: generate random initials
        % 随机状态个数
        num = 1000; %num of random initial states
        % 噪音强度
        alpha = 0.01; %noise strength in random initial states

        % 生成一堆起始细胞表达状态       
        randomXInits=generateRandomInitialStates(num,alpha,hopland);        
        
        %生成每个基因表达水平随细胞状态的变化矩阵 和 相应的权重矩阵
        % Step 3: generate test traj from real data
        [realTraj,weight]=generateTraj(hopland);        
        
        %参数优化
        % Step 4: parameter optimization
        maxIts=2000;

        paramInit=parameterOptimization(maxIts,randomXInits,fittingData,fittingDataTemp,hopland,realTraj,weight);    
        % 将其赋值给paramInit    
        hopland.paramInit=paramInit;
        
    %% clustering
    else
        %如果没有时间序列信息，就直接用初始化参数
        paramInit=initializeParam(hopland.orgiFitData);
    end


    %地形构建
    %% Step 5: landscape construction
    % 得到GPLVM的模型，1*438的细胞能量，31*31所有网格点的能量矩阵，两个对应的二维网格坐标length(y)*length(x)
    [model,energyLand,ENERGYLAND,X,Y]=constructLandscape(paramInit,fittingData,hopland);
    % 赋值进这个结构体
    hopland.model=model;
    hopland.ENERGYLAND=ENERGYLAND;
    hopland.energyLand=energyLand;
    hopland.X=X;
    hopland.Y=Y;

  
    %伪时估计
    %% Step 6: pseudotime estimation
    display=1;
    % 输出（根到各个细胞的距离1*438，不同细胞所属阶段 和 根到不同细胞距离 之间的相关系数1*1）
    [dist,coef]=calculateDistance(hopland,display,givenStartPoint,ifdoComparison);
    % 增加结构体参数
    hopland.dist=dist;
    hopland.coef=coef;
end