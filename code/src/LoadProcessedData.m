%函数LoadProcessedData
%参数（数据集的名字）
%返回一个结构体hopland
%变量（开始阶段的细胞index列表，每个阶段的细胞数，选择的基因名字，选择基因的表达矩阵，发育时间轴，所有细胞的标签，数据集本身）

function hopland=LoadProcessedData(dataset)
    %判断dataset的名字
    if( strcmp(dataset,'guo2010') )
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessRobinson(1,0);
        %通过观察，设置为1的原因是细胞标签比较细化，和时间比较相关，这样的话在挑选基因的时候就用t检验
        hopland.ifTimeseries=1;

    elseif( strcmp(dataset,'deng2014_1') )
        %include fibroblast and adult liver （成纤维细胞和成年人肝脏细胞） 
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessDeng(1,1);
        hopland.ifTimeseries=1;

    elseif( strcmp(dataset,'yan2013') )
        %exclude fibroblast and adult liver
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessYan(1,1);
        hopland.ifTimeseries=1;        

    elseif( strcmp(dataset,'SyntheticData') )
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessSyn(0,0);
        hopland.ifTimeseries=1;

    %从这里开始的数据并没有时间的标签，这样分类就比较粗，选取基因的时候直接用var来选  
    elseif( strcmp(dataset,'LPS') )
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessLPS(0,1);
        hopland.ifTimeseries=0;
        
    elseif( strcmp(dataset,'HSMM') )
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessHSMM(0,1);
        hopland.ifTimeseries=0;
        
    elseif( strcmp(dataset,'ES_MEF') )
        [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessESMEF(0,0);
        hopland.ifTimeseries=0;
    end
    
    
    hopland.startRefRange=startRefRange;
    hopland.cellStates=cellStates;
    hopland.selectedGeneNames=selectedGeneNames;
    hopland.orgiFitData=filteredExpData;
    hopland.developLine=developLine;
    hopland.cellLabels=cellLabels;
    hopland.dataset=dataset;

end


%ProcessHSMM;
% ProcessDeng;
% ProcessSimulate;