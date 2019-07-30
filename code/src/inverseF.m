% x is a float
% 函数inverseF
% 输入（细胞的某个基因表达值1*1，整体基因的均值和标准差2*48，基因序号）
% 输出（1*1的inverse结果） 
function invfx=inverseF(xValue,fittingData,genei)
    %old
    invfx = -fittingData(2,genei)*log(1./xValue-1)+fittingData(1,genei); %1/(1+exp(-(x(gene)-mixedmu)/mixedtheta));
    %new
%     syms x
%     invfx=double(subs(fittingDataInverse{genei},x,xValue));
end