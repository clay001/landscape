function [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessLPS(ifNormaliseData,ifFilterGenes,SelectedGenesIndex)
    %load data 
    load('./sample_data/lps/lps.mat');

    % genes, labels, geneExpData, cellStates
    count=tabulate(pro.cell_stage);
    %统计出每个阶段的频数
    cellStates=count(:,2);
    %取出基因表达矩阵
    geneExpData=pro.expr;
    %开始阶段的基因下标
    startRefRange=1:cellStates(1);

    numGenes=length(pro.gnames);
    developLine=0;
    cellLabels=pro.cell_stage;
    
    %% Step 0: Process data
    %同样的归一化操作
    if ifNormaliseData
        geneExpData = zscore(geneExpData);
    end
    
    %在过滤基因这里操作不一样了
    %filter genes
    if ifFilterGenes
        %基因表达数据每一列的方差降序排列，I对应的是该列在原矩阵中的位置
        [Y,I]=sort(var(geneExpData),'descend');
        %选择top的基因
        %top选择的原则是最大差异的基因在原排列中的位次除以5，感觉有点缺乏说服力，不过也行
        SelectedGenesIndex=I(1:floor(I/5));
    end
    
    if ~exist('SelectedGenesIndex')
        selectedIndex=1:numGenes;
    else
        selectedIndex = SelectedGenesIndex;
    end
    
    filteredExpData=geneExpData(:,selectedIndex);
    selectedGeneNames=pro.gnames(selectedIndex);

end
