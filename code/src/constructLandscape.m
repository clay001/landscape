% 函数constructLandscape
% 输入（模型参数，包含每个基因的均值和标准差的double数组2*48，包含数据的结构体）
% 输出（GPLVM的模型，1*438的细胞能量，31*31所有网格点的能量矩阵，两个对应的二维网格坐标length(y)*length(x)）
function [model,energyLand,ENERGYLAND,X,Y]=constructLandscape(paramInit,fittingData,hopland)
    orgiFitData=hopland.orgiFitData;
    
    %learning model
    % 返回模型，原始表达矩阵
    [model,Y]=runGPLVM(orgiFitData,'linear','isomap',20); %isomap

    %% plot mapping
    % 关掉画模型的两个窗口
    close all;
    %plotMappingResult(model,cellStates,1,1); %style=1
%     hopland.model=model;
%     plotMappingResult(hopland,1,1);

    %% land energy
    % m为48
    m=size(orgiFitData,2);
    % 权重矩阵先设为全零
    weightMatrix= zeros(m);
    count=1;
    for j=1:m
        for i=1:j
            weightMatrix(i,j)=paramInit(count);
            weightMatrix(j,i)=weightMatrix(i,j);
            count=count+1;
        end
    end
    % 依次取出填充

    % 再依次取出A，sigma，W(权重矩阵改名)
    A=paramInit( ((m*m+m)/2+m+1):((m*m+m)/2+2*m));
    sigma=paramInit( ((m*m+m)/2+m*2+1):((m*m+m)/2+3*m))';   
    W=weightMatrix;
    
    % 如果时间标签不足，sigma设0 （不是本来就是0吗？）
    if hopland.ifTimeseries == 0 
        sigma=0;
    end
    % sigma横着放（因为取出来的时候由于它是一维的，默认都会当成向量竖着放）
    if size(sigma,1) ~= 1
        sigma=sigma';
    end

    % model.X是438*2 , 应该是438个细胞降维之后的两个component
    % YpredReal是438*48
    YpredReal = gpOut(model, model.X);
    energyLand=zeros(1,size(YpredReal,1));
    % ii从1到438
    for ii=1:size(YpredReal,1)
        temp=YpredReal(ii,:);
        % 将每个细胞的表达值做非线性变换， 1*48
        fx=F(temp,fittingData);
        
        invFx=0;
        % 并行for循环
        parfor genei=1:size(YpredReal,2)
            % 对每一个基因
            % 先定义一个句柄，把这三个参数传到inverse里
            f = @(xGene,fittingData,genei)inverseF(xGene,fittingData,genei);
            % integral函数fun，xmin,xmax
            % 这一部分算积分项，第一项定义一个函数，里面的x每次由后面的fx喂
            invFx = invFx + A(genei)*integral(@(x)f(x,fittingData,genei), 0, fx(genei)); %invfx=inverseF(x,fittingData);
        end
        % 计算能量，填入1*438的energyLand里
        energyLand(ii) = -0.5 * fx * W * fx' + invFx - sum(fx.*sigma);
    end
    

    %% generate mesh landscape
    % 从模型的两个component里取出大小值
    topMax=max(model.X);
    minMin=min(model.X);
    % 
    n=30;
    i_start = minMin(1)- (topMax(1)-minMin(1))/4;
    i_end= topMax(1) + (topMax(1)-minMin(1))/4;

    j_start= minMin(2) - (topMax(2)-minMin(2))/4;
    j_end = topMax(2) + (topMax(2)-minMin(2))/4;

    % 画二维网格，[X,Y] = meshgrid（x,y）X是一个矩阵，每一行是x的一个副本
    % Y也是一个矩阵，每一列是y的一个副本
    % 返回一个二维网格坐标
    [X,Y] = meshgrid(i_start:(i_end-i_start)/n:i_end, j_start:(j_end-j_start)/n:j_end);
    
    % Xpred大小为961*2 
    Xpred=zeros((n+1)*(n+1),2);k=1;
    % 分别等分出30个子块，实际上有31对值
    for i=i_start:(i_end-i_start)/n:i_end
        for j=j_start:(j_end-j_start)/n:j_end
            % 依次填入i和j
            Xpred(k,:)=[i,j];
            k=k+1;
        end
    end
    

    % Ypred得到 961*48
    Ypred = gpOut(model, Xpred);
    % ENERGYLAND大小设为31*31
    ENERGYLAND=zeros(n+1,n+1);

    for ii=1:n+1
        parfor jj=1:n+1
            % temp取出的是网格中的index号，逐列编号
            % 取出这个号，也就是这个点的表达情况1*48
            temp=Ypred((jj-1)*(n+1)+ii,:);
            % 做一个非线性变换
            fx=F(temp,fittingData);
            invFx=0;
            % 相当于是对每个点都计算一个能量
            for genei=1:size(Ypred,2)
                f = @(xGene,fittingData,genei)inverseF(xGene,fittingData,genei);
                invFx = invFx+A(genei)*integral(@(x)f(x,fittingData,genei), 0, fx(genei)); %invfx=inverseF(x,fittingData);
            end           
            ENERGYLAND(ii,jj) = -0.5 * fx * W * fx' + invFx - sum(fx.*sigma);
        end
    end

    %最后得到的大写energy应该是31*31的矩阵，包含网格所有点的能量
    
end