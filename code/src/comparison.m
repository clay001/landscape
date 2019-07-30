% 函数comparison
% 输入（数据结构体）
% 输出（打印两个值，返回不同细胞所属阶段 和 根到不同细胞距离 之间的相关系数）
function coef = comparison(hopland)
    % 值传递
    reference=hopland.dataset;
    dist=hopland.dist;
    cellLabels=hopland.cellLabels;

    % 判断数据集的名字
    if( strcmp(reference,'guo2010') )
        % reference重新定义为1*438的double，每个值对应的是每个细胞处在的细胞阶段
        reference=[ones(1,9),2*ones(1,19),3*ones(1,23),4*ones(1,44),5*ones(1,75),6*ones(1,109),7*ones(1,159)];
    elseif( strcmp(reference,'deng2014_1') )
        reference=[ones(1,4),1.7*ones(1,8),2*ones(1,20),2.5*ones(1,10),...
            3*ones(1,14),4*ones(1,47),...
            5*ones(1,58),6*ones(1,43),6.5*ones(1,60),7*ones(1,30),8*ones(1,23)];      
    elseif( strcmp(reference,'deng2014_2') )
        reference=[ones(1,4),1.7*ones(1,8),2*ones(1,20),2.5*ones(1,10),...
            3*ones(1,14),4*ones(1,47),...
            5*ones(1,58),6*ones(1,43),6.5*ones(1,60),7*ones(1,30)];   
    elseif( strcmp(reference,'yan2013') )
        reference=[ones(1,6),2*ones(1,6),3*ones(1,12),4*ones(1,20),5*ones(1,16),6*ones(1,30),7*ones(1,34)];
    elseif( strcmp(reference,'SyntheticData') )
        reference=cellLabels;
    elseif( strcmp(reference,'LPS') )
        reference=cellLabels;
    elseif( strcmp(reference,'HSMM') )
        reference=cellLabels;
    elseif( strcmp(reference,'ES_MEF') )
        reference=cellLabels;
    end
    
    %% comparison
    % D是1*438，表示根到不同细胞的距离
    D=dist;
    % 找到dist是负的位置，一般是没有的
    index=find(dist==-1);
    % setdiff（A,B）函数求两个数组的差集
    % 返回A中存在但B中不存在的数据，不包含重复项，是有序的
    % reference是1*438
    reference = reference(setdiff(1:length(reference),index));
    % D取出D中不为距离负数的部分, 1*438
    D = D(setdiff(1:length(D),index));
    % corrcoef函数求线性相关系数，放进两列求系数 2*2 （原矩阵i列和j列的相关系数，应该是对称的）
    temp=corrcoef(reference',D');
    % coef取出（1，2）的元素
    coef=temp(1,2);
    
    % lgd为负数位置的字符串，一般为空
    lgd = num2str(index,'#%d ');
    % 打印出相关系数和不能连接的点
    details1=strcat('The correlation coefficients: ',num2str(coef));
    details2=strcat('Points that cannot be connected: ',lgd);
    sprintf('%s\n%s\n',details1,details2)
end
