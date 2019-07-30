% 函数hopfieldNetworkContinuousModel
% 输入（这里的权重指的是根据相关系数填入的对称阵，tTime，sigma, A,sigmaW，初始细胞表达水平，每个基因的均值和标准差，发育线
% 输出（按发育尺取出了模拟的Xem从初始状态开始走过所有细胞阶段的表达情况
function [Xem]=hopfieldNetworkContinuousModel(W, T, sigma, A,sigmaW, Xzero, fittingData, developLine)
    %% dynamics
    % 正态分布的随机数，建议改用rng(100,'v5normal')
    randn('state',100);
    %M = size(W,1);

    % 对称W点乘上sigmaW
    W=W.*sigmaW;

    %% SDE
    % N = T*64; dt = T/N;
    % dW = sqrt(dt)*normrnd(0,1,M,N); % Brownian increments
    %
    % R = 4; L = 442; %N/R; % L EM steps of size Dt = R*dt
    % Xem = zeros(L+1,M); % preallocate for efficiency
    % Xem(1,:)=Xzero;
    % Xtemp = Xzero;
    % for j = 1:L
    %     for gene=1:M
    %         Winc = sum(dW(gene, R*(j-1)+1:R*j));
    %         Xtemp(gene) = Xtemp(gene) -A(gene)*Xtemp(gene) + W(gene,:)*F(Xtemp,fittingData)' + sigma*Winc;
    %         Xem(j+1,gene) = Xtemp(gene);
    %     end
    % end

    %% ODE
    % ff = @(t,Xtemp) tempF(t,Xtemp,fittingData,M,W,sigma,A);
    % [t,Xem]=ode45(ff,0:T/100:T,Xzero);      %
    % plot(t,Xem)


    %% SSA
    % [t,Xem] = ssa_example(Xzero,fittingData,W,sigma,A,T);


    %% Euler method
    % 这里T起到控制h的作用，干什么的不是很清楚
    h=0.1;
    % 按1传进来的话就是10<50，h直接改为0.02
    if T/h<50
        h=T/50;
    end
    % 估计状态Xem，这样就开了51个状态
    Xem=zeros(T/h+1,size(W,2));
    % 第一个直接还是原本的细胞状态
    Xem(1,:)=Xzero;
    count=2;
    % 循环估计
    for i=2:T/h+1
        % 在上一个的基础上加上 h和变化率的乘积，当作下一个Xem的状态（欧拉方法）
        Xem(count,:) = Xem(count-1,:) + h*tempF([],Xem(count-1,:)',W,sigma,A,fittingData)';
        count=count+1;
    end

    % 把发育尺扩大预测状态的细胞个数倍  再向下取整
    % 这里得到的等分为细胞阶段数目的index抽样
    newindex=floor(developLine*size(Xem,1));       %[1:count/(size(realTraj,2)-1):size(Xem,1),size(Xem,1)];
    % 再人为地将第一个元素设为1，因为matlab的角标是从1开始的
    newindex(1)=1;
    % 这里再floor一下好像有点多余，就是按发育尺取出了Xem等间距细胞表达情况
    Xem=Xem(floor(newindex),:);
end





% 底下定义一个函数tempF
% 传入t是空的，Xtemp为上一个细胞状态的转置，对称权重阵，修正项sigma,衰减项A,每个基因的均值和标准差
% 传出一个f，48*1的变化率
function f=tempF(t,Xtemp,W,sigma,A,fittingData)
    % F函数写在另一个文件中，喂入上一细胞状态的转置48*1和每个基因的均值和标准差2*48
    % 返回的是 1*48的激活函数作用后的加权和
    % result又转置，成为48*1
    resultF = F(Xtemp,fittingData)';
    % 把A中较大的维度躺着放
    if size(A,1)>size(A,2)
       A=A'; 
    end
    
    % 如果resultF躺着就把他竖过来
    if size(resultF,1)==1
        resultF=resultF';
    end
    
    % sigma横过来
    if size(sigma,1)>1
        sigma=sigma';
    end
    
    % 总之就是为了这个公式好乘
    f=-A'.*Xtemp + W*resultF + sigma';  
end

