% 函数generateRandomInitialStates
% 输入（生成点个数，噪音强度，原始数据结构）
% 输出（随机抽取并且加上噪音扰动的细胞及其表达水平）
function randomXInits=generateRandomInitialStates(num,alpha,hopland)
    orgiFitData=hopland.orgiFitData;
    startRefRange=hopland.startRefRange;
    
    % mean mu ???
    % 求起始细胞状态的所有基因表达数据的标准差
    mu=std(orgiFitData(startRefRange,:));

    randomXInits = zeros(num,size(orgiFitData,2));
    % 开始循环生成
    for i=1:num
        % randi均匀分布的伪随机整数,个数是起始细胞状态的数量
        origNum = randi(length(startRefRange));
        % 取出随机生成的细胞
        Xzero = orgiFitData(origNum,:); %randn(1,size(leftgeneExpData,2)); % initial values

        % 进行填充，在随机抽样的细胞数据值上
        % rand（x）均匀分布的随机整数
        % 前半部分是直接加上一个小的噪音
        % 后半部分alpha和mu前面干的事是随机生成一个值为正负1的矩阵进行扰动
        randomXInits(i,:) =Xzero+0.001*rand(1,size(Xzero,2))+((rand(1,size(orgiFitData,2))>0.5)*2-1)*alpha.*mu;
    end

end