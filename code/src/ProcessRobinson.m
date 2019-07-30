%函数ProcessRobinson
%参数(是否归一化，是否自动选取基因，预指定基因index)
%返回（开始阶段的细胞index列表，每个阶段的细胞数，选择的基因名字，选择基因的表达矩阵，发育时间轴，所有细胞的标签）
%猜测标签的命名意义是多少细胞阶段，以及几月几号取出的时间
function [startRefRange,cellStates,selectedGeneNames,filteredExpData,developLine,cellLabels]=ProcessRobinson(ifNormaliseData,ifFilterGenes,SelectedGenesIndex)
    %load data 
    load('./sample_data/guo2010/guo2010PData2.mat');

    % genes, labels, geneExpData, cellStates

    % tabulate函数的功能是创建向量信息数据频率表，第一列是参数，第二列是频数，第三列是频率
    count=tabulate(pro.cell_stage);

    %用cellStates取出第二列，也就是频数信息，参数就可以隐含在index里
    cellStates=count(:,2);

    %geneExpData取出基因表达的矩阵
    geneExpData=pro.expr;

    %返回的第一个参数是第一个细胞阶段的数目
    startRefRange=1:cellStates(1);

    %基因的数量，robinsob的话就是48，是一个1*1的double
    numGenes=length(pro.gname);

    %起止和python不太一样，是一个双头闭的区间
    %developLine把所有的cellStates都等分在0-1区间内
    developLine=0:  1/(length(cellStates)-1)  :1;

    %cellLabels把数据中的cell特征直接取出来
    cellLabels=pro.cell;
    

    %% Step 0: Process data
    %如果需要归一化数据
    if ifNormaliseData
        %调用了zscore函数，做的事情就是减去均值之后再除以标准差，常用的操作
        geneExpData = zscore(geneExpData);
    end
    

    %filter genes
    %如果需要过滤基因
    if ifFilterGenes
        %调用cunsum函数进行元素累加，默认后面参数置为1（按列累加），返回的是一个和原来形状相同的矩阵，每行表示累加到当前行数
        %时候的值
        %结果进行转置存储在ends里，表达的信息到第几个阶段为止的所有细胞总数，在程序中表述成每个阶段结束的细胞index
        ends=cumsum(cellStates)';
        
        %在matlab里可以直接利用end进行切片控制，取出前六个细胞阶段，并且加上一个初始状态细胞index(设为1)
        %在程序中表述为每个阶段开始的细胞index
        starts=[1,ends(1:end-1)+1];

        %基因表达矩阵的大小是438*48，取出的是一个1到48的索引表
        selectedIndex=1:size(geneExpData,2);


        %从1到6的细胞阶段循环
        for i=1:length(cellStates)-1
            %里循环是从外循环的后一个细胞阶段开始
            for j=i+1:length(cellStates)
                %调用了mattest函数，格式为[PValues, TScores, DFs] = mattest(DataX, DataY)
                %输入对比的是不同阶段的细胞基因表达矩阵（为了行数匹配，进行转置，用基因名字作为行）
                [pvalues, tscores] = mattest(geneExpData(starts(i):ends(i),:)', geneExpData(starts(j):ends(j),:)');
                %得到fdr,positive false discovery rate和qvalues(hypothesis testing error)
                [pFDR, qvalues] = mafdr(pvalues);
                %intersect返回参数的共有数据（排序且不包含重复项），find函数是返回一个包含参数中非零元素索引的向量
                %在索引表中取出qvalue值小于阈值的
                selectedIndex=intersect(selectedIndex,find(qvalues<0.05));
            end
        end

        %上述相当于做了一个stage间的两两test，根据结果不断地缩小索引表的范围，最后得到选中的index
        %也就是差异性最大的基因们

    %如果不要求程序来指定index，又没指定gene的index，就默认为全体gene
    elseif ~exist('SelectedGenesIndex')
        selectedIndex=1:numGenes;
    %如果指定了那就直接用指定的index
    else
        selectedIndex = SelectedGenesIndex;
    end
    %可以看到如果需要程序来指定index，又给了指定index的话，指定的要求是会被忽略的

    %返回过滤后的所有细胞的特定基因的表达矩阵
    filteredExpData=geneExpData(:,selectedIndex);
    %同时也返回特定基因的名字列表，也就是上述矩阵的列名
    selectedGeneNames=pro.gname(selectedIndex);

end
