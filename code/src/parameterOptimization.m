% 函数parameterOptimization
% 输入（循环次数，随机起始数据，每个基因的均值和标准差，fit高斯模型的参数，包含数据的结构体，基因轨迹和标准差权重48*7）
% 输出（经过迭代之后的最优参数1*1320）
function paramInit=parameterOptimization(maxIts,randomXInits,fittingData,fittingDataTemp,hopland,realTraj,weight)
    %显示变量的值（和直接打印没差别）
    disp('Parameter optimization...');

    orgiFitData=hopland.orgiFitData;
    developLine=hopland.developLine;
    
    %initilize parameters
    %初始化参数
    paramInit=initializeParam(orgiFitData);

    %% adaption
    % 学习率
    learnRate = 0.3;
    % 基因个数
    N=size(orgiFitData,2);
    % 真实基因矩阵
    xReal=orgiFitData;
    % 最终结果
    savedResults=zeros(1,2);
    % 最终参数
    savedParam=zeros(1,length(paramInit));
    % 存储index
    cc=1;

    % 得到总误差，和模拟的结果
    [g1,xSimulate1,xSimulate,xSimulate2] = trySN(paramInit,randomXInits,fittingData,hopland,realTraj,weight);
    % 简单的赋值
    g0=g1;

    % 开始迭代
    for its=1:maxIts    
        % 不断地从起始值模拟，g1得到trajDiff
        [g1,xSimulate1,xSimulate,xSimulate2] = trySN(paramInit,randomXInits,fittingData,hopland,realTraj,weight);
        % 基于模拟值和高斯模型计算cdfDiff，我很奇怪这里传入的为什么是simulate1
        x2=mycdfDiff(xSimulate1,xReal,realTraj,fittingDataTemp);

        % 把这两个误差传入savedResults
        savedResults(cc,:)=[g1,sum(x2)];
        % 先传入初始参数
        savedParam(cc,:)=paramInit;
        % 挪到下一个index，这里可以直接通过增加cc来向之前的save里面添加行
        cc=cc+1;


        %adapt learning rate
        % 如果traj误差变大了
        if g1>g0
            % 降低学习率
            learnRate=0.5*learnRate;
            paramInit = paramInitbak1;
            % 直接下一轮迭代
            continue;
        % 如果误差变小了
        elseif g1<g0
            % 加大学习率
            learnRate=(1+0.02)*learnRate;
            % 结果满意，继续运行
        end

        % 如果两次误差非常相近，且两个误差又不相同
        if g0-g1<10^(-5) && g0~=g1
            %init设为back
            paramInit=paramInitbak1;
            % 并且结束整个循环
            break;
        end
        % g0表征上一轮的error，g1表征当前的error
        g0=g1;

        % 记录上一回合的参数
        paramInitbak1 = paramInit;

        % 如果学习率太小
        if learnRate<10^(-6)
            % 回退上一次的参数
            paramInit=paramInitbak1;
            % 并且结束整个循环
            break;   
        end
        
        % 48*48的权重矩阵      
        weightMatrix= zeros(size(xReal,2),size(xReal,2));

        count = 1;
        % 取出参数对矩阵进行填充
        for j=1:size(xReal,2)
            for i=1:j
                weightMatrix(i,j)=paramInit(count);
                weightMatrix(j,i)=weightMatrix(i,j);
                count=count+1;
            end
        end

        % 取出参数sigmaW
        sigmaW =paramInit(((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+1):((size(xReal,2)*size(xReal,2)+size(xReal,2))/2+size(xReal,2))); %
        Ttime=1;

        % 类似hopfilednetwork里的操作
        % % 按1传进来的话就是10<50，h直接改为0.02
        h=0.1;
        if Ttime/h<50
            h=Ttime/50;
        end
        
        % 细胞阶段个数
        numPoints=length(developLine);    
        newTime=1:numPoints;

        % allF初始化为7*48的全0矩阵
        allF=zeros(size(xSimulate2));

        % 对每个阶段
        for i=1:size(xSimulate2,1)
            % F函数，喂入1*48 和 2*48 
            % 得到 1*48
            allF(i,:)=F(xSimulate2(i,:),fittingData);
        end
        % 最后allF是7*48
        
        % 1*48
        da_i=zeros(1,N);        
        dI_i=zeros(1,N);
        dC_i=zeros(1,N);

        dw_ij=zeros(1,(N*N+N)/2);

        count = 1;  
        for genei=1:N
            % 取出7行矩阵的表达值7*1 
            % 转置变成1*7
            temp=xSimulate2( newTime,genei)';

            % 梯度下降法对A，I，C求导
            ai=0;ii=0;ci=0;
            for i=1:numPoints
                if i==1
                    iii=2;
                else 
                    iii=i;
                end
               ai = ai + weight(i)* (realTraj(genei,i)-temp(i)) * h * sum(xSimulate2(1:(newTime(iii)-1),genei));                          
               ii = ii + weight(i)* (realTraj(genei,i)-temp(i)) * h * (newTime(iii)-1);
               if newTime(iii) - 1 == 1
                    ci = ci + weight(i)* (realTraj(genei,i)-temp(i)) * h * (weightMatrix(genei,:)*(allF(1:(newTime(iii)-1),:))');
               else
                   ci = ci + weight(i)* (realTraj(genei,i)-temp(i)) * h * (weightMatrix(genei,:)*sum(allF(1:(newTime(iii)-1),:))');
               end  
            end
            da_i(genei)=-2/N*ai;            
            dI_i(genei)=2/N*ii;
            dC_i(genei)=2/N*ci;


            % 对W求导
            for j = 1:genei
                wij = 0;
                for i=1:numPoints
                    if i == 1
                        iii=2;
                    else
                        iii=i;
                    end
                    wij = wij + weight(i) * (realTraj(genei,i)-temp(i)) * h * sigmaW(genei) * sum(allF(1:(newTime(iii)-1),j));
                end
                % 填入
                dw_ij(count)=2/ N * wij;
                count=count+1;
            end 
        end

        % 按照梯度下降更新参数1320*1
        paramInit1 = paramInit+learnRate*[dw_ij,dC_i,da_i,dI_i]';
        k=(N*N+N)/2;

        % find w之后，这些值里小于0的值的index
        negIndex=find(paramInit1((k+1):end)<=0)+k;
        % 这些小于0的位置还是换回原来的参数值
        paramInit1(negIndex)=paramInit(negIndex);
        % 把这一轮的参数值传递给上一轮的参数值，进行下一轮迭代
        paramInit=paramInit1;

    end

    %% 选取的时候选取使得误差和最小的参数值
    [a,I]=min(sum(savedResults'));
    % 返回这个序号对应的参数值
    paramInit=savedParam(I,:);
    disp('done');

end
