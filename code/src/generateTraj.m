% 函数generateTraj
% 输入（包含数据的结构体）
% 输出（命名为realTraj的每个基因随细胞状态变化的表达均值48*7  和  每个基因在每个细胞阶段的权重，表征表达水平的统一程度48*7）
function [realTraj,weight]=generateTraj(hopland)

    orgiFitData=hopland.orgiFitData;
    cellStates=hopland.cellStates;
    % 时间点的个数设为细胞状态的个数
    numTimePoints=length(cellStates);

    % 真实的轨迹大小为基因个数*细胞状态个数
    realTraj=zeros(size(orgiFitData,2),numTimePoints);
    % 权重同形
    weight=zeros(size(orgiFitData,2),numTimePoints);

    %开始循环
    for i=1:size(orgiFitData,2)
        % 对每一个基因
        geneID=i;
        % 对第一个初始阶段，这里用cumsum不知道有什么意义，总之就是对每个基因的第一阶段的细胞求均值
        % 填充进该行第一个位置
        realTraj(i,1)=mean(orgiFitData(1:cumsum(cellStates(1)),geneID));
        % 权重是标准差的倒数，也就是说数据之间差异性越小，权重越大
        weight(i,1)=1./std(orgiFitData(1:cumsum(cellStates(1)),geneID));

        for j=2:numTimePoints
            % 中间这一堆做的事情就是取出数据中的相对应阶段的细胞index
            % 然后还是一样填充均值和标准差
            realTraj(i,j)=mean(orgiFitData((sum(cellStates(1:(j-1)))+1):sum(cellStates(1:j)),geneID));
            weight(i,j)=1./std(orgiFitData((sum(cellStates(1:(j-1)))+1):sum(cellStates(1:j)),geneID));
        end
    end

    % 控制权重上界
    weight(weight>10)=10;
end