% 函数initializeParam
% 输入（表达数据）
% 输出（初始化的所有参数1320*1）（包括基因两两之间的相关系数和sigmaW,a,sigma）
function paramInit=initializeParam(orgiFitData)
    %initilize parameters
    % 取出基因的个数
    tempn = size(orgiFitData,2);
    % 参数的总个数
    params = zeros((tempn^2+tempn)/2+tempn*3,1); %

    % corr可以比较两个矩阵之间得到线性和秩相关性，也可以对单个矩阵返回列之间两两的相关系数
    % 应该是48*48
    W=corr(orgiFitData);
    % triu取出矩阵的上三角部分,不等于0的全部设为1
    % 再用W取出来的数值是逐列从上三角矩阵不为0的元素开始取出来
    % 其实这里reshape默认了不会出现等于0的值
    temp=reshape(W(triu(W)~=0),1,(tempn^2+tempn)/2);

    % 对于相关性矩阵中的每一个参数
    for i=1:(tempn^2+tempn)/2
        % 依次填入参数
        params(i)=temp(i);
    end

    %从这之后还增开3*tempn的参数，全部填入1，1，0
    %sigmaW
    for i=((tempn^2+tempn)/2+1):((tempn^2+tempn)/2+tempn)
        params(i)=1; %abs(randn(1)); 
    end
    %a
    for i=((tempn^2+tempn)/2+tempn+1):((tempn^2+tempn)/2+tempn*2)
        params(i)=1; %abs(randn(1));
    end
    %sigma
    for i=((tempn^2+tempn)/2+tempn*2+1):((tempn^2+tempn)/2+tempn*3)
        params(i)=0; %abs(randn(1));
    end

    paramInit=params;

end