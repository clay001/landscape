%函数fitMixtureGaussian
%输入（一个包含数据的结构体）
%返回（包含每个基因的均值和标准差的double数组2*48和144*4的数组纵列是component,横列三个一组是模型的参数，mu和sigma）
function [fittingData,fittingDataTemp]=fitMixtureGaussian(hopland)
    orgiFitData=hopland.orgiFitData;
    
    %返回的数据先重设一个全零的double数组
    fittingData=zeros(2,size(orgiFitData,2));
    %第一行用原始数据的平均值填充
    fittingData(1,:)=mean(orgiFitData);
    %第二行用原始数据的标准差填充
    fittingData(2,:)=std(orgiFitData);
    
    %设一个大小为1*gene个数的空元胞
    fittingDataTemp=cell(1,size(orgiFitData,2));
    
    %随机种子
    seed = 1;
    %创建随机数流
    s = RandStream('mt19937ar','Seed',seed);
    %创建一个统计结构
    opts = statset('Streams',s, 'MaxIter', 2000);

    %对每一个基因
    for geneID=1:size(orgiFitData,2)        
        %aic和bic先置为全零的数组
        aicScores = zeros(1,4);
        bicScores = zeros(1,4);
        %reslut是一个元胞，用于存储高斯模型
        results=cell(1,4);

        for k=1:4
            %fitdata用某一列基因表达水平， opts，非负的数字加在协方差矩阵的对角使其正定
            %现在大多都用fitgmdist(orgiFitData(1:end,geneID),k)       
            gmm = gmdistribution.fit(orgiFitData(1:end,geneID),k,'Options',opts,'Regularize',0.000001);
            %分别对外部变量进行填充
            aicScores(k) = gmm.AIC;
            bicScores(k) = gmm.BIC;
            results{k}=gmm;
        end
        %返回最小的AIC和它的k
        [a,k]=min(aicScores);
        
        %check
        %设置一个阈值
        threshold=0.05;
        %取出该模型
        gmm=results{k};

        %看看这时候的每个参数值是否存在小于阈值的情况，如果小于阈值了，就尝试减小k的值
        if k==2
            if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold
                k=1;
            end
        elseif k==3
            if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold || gmm.PComponents(3)<threshold
                k=2;
                gmm=results{2};
                if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold
                    k=1;
                end
            end
        elseif k==4
            if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold || gmm.PComponents(3)<threshold || gmm.PComponents(4)<threshold 
                k=3;
                gmm=results{3};
                if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold || gmm.PComponents(3)<threshold
                    k=2;
                    gmm=results{2};
                    if gmm.PComponents(1)<threshold || gmm.PComponents(2)<threshold
                        k=1;
                    end
                end
            end
        end
        
        %final fitting
        
        %在该ID的元胞内填入相应的参数、均值、协方差，不满4个用0填充
        if k<4
            %这里有点confuse,以我调试的情况为例，k是2，首先第一行放component,剩余用0填充
            %第二行放mu（1*k）的转置，剩余用0填充
            %第三行展开sigma（1*1*k），剩余用0填充
            %zeros是因为这里只fit了一列数据，所以前面的维数都是1*1...*1
            fittingDataTemp{geneID}=[results{k}.PComponents,zeros(1,4-k);
                results{k}.mu',zeros(1,4-k);
                results{k}.Sigma(:,:),zeros(1,4-k)];
        else        
            fittingDataTemp{geneID}=[results{k}.PComponents;
                results{k}.mu';
                results{k}.Sigma(:,:)];
        end       
        
        %打印基因名字和k值，显示正在处理
        fprintf('%s Gene#%i: k=%i\n','Processing ',geneID,k);
    end
    
    %将元胞类型转化为基础的数组,程序运行完以后1*48的元胞内都是一个3*4的double
    %相同列数拼接成一个3*48=144 144*4的double
    fittingDataTemp=cell2mat(fittingDataTemp');

end