% 函数calculateDistance
% 输入（数据结构体，是否显示，给定初始点（1），是否比较（1））
% 输出（根到各个细胞的距离1*438，不同细胞所属阶段 和 根到不同细胞距离 之间的相关系数1*1）
function [dist,coef]=calculateDistance(hopland,display,StartPoint,ifdoComparison)
    % 参数传递
    model=hopland.model;
    cellStates=hopland.cellStates;
    energyLand=hopland.energyLand;

    % x和y都是438*1
    % 表征438个细胞的两个component，如果设置再原点上的等比例值
    x=(model.X(:,1)-min(model.X(:,1)))/(max(model.X(:,1))-min(model.X(:,1)));
    y=(model.X(:,2)-min(model.X(:,2)))/(max(model.X(:,2))-min(model.X(:,2)));
    % z为1*438, 也是等比例放缩
    z=2*(energyLand-min(energyLand))/(max(energyLand)-min(energyLand));    
    
    % manifold保存这三个值， 438*3
    manifold=[x,y,z']; 
    
    %% construct mesh
    % 三维上进行三角剖分算法
    % delaunayTriangulation函数，输出的是一个delaunayTriangulation的数据类型
    % 包含约束边dt.Constraints（空），三角剖分点dt.Points（438*3）和三角剖分连接列表dt.ConnectivityList（2561*4）
    dt = delaunayTriangulation(manifold);
    % 这个就和manifold里的值是一样的，转置后为3*438，表征每个点的坐标
    vertex=dt.Points';
    
    % dt1用原始未经处理的坐标进行剖分 2478*4
    dt1 = delaunayTriangulation([model.X(:,1),model.X(:,2),energyLand']);
    % 3*438
    vertex_orgi=dt1.Points';

    % 二维的三角剖分859*3
    dt = delaunayTriangulation([x,y]);
    % 面，连接列表为859*3
    % 每行代表一个三角形，元素是顶点ID
    faces=dt.ConnectivityList;

    %% set start point
    % 如果不存在初始点
    if ~exist('StartPoint')
        % 那就让初始点为细胞状态中最大的，a是最大的值，StartPoint为它的细胞序号
        [a,StartPoint]=max(energyLand);
    end

    %% cal dist
    % distanceMatrix大小为438*438
    distanceMatrix=zeros(size(vertex,2),size(vertex,2));
    % 对每一个细胞
    for start_points=1:size(vertex,2)
        % 设定当前初始点位置
        options.start_points=start_points;
        % 结束点位置设空
        options.end_points = [];
        
        % [D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options)
        % D 返回该细胞到其他所有细胞的距离
        % 三个返回都是438*1
        [D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points);
        % 填入distanceMatrix
        distanceMatrix(start_points,:)=D;
    end
    
    % 最小生成树图
    % 函数graphminspantree，输入距离矩阵，给定root，返回一颗稀疏的最小生成数和每个点的前任点
    % 这里距离矩阵在输入的时候进行sparse处理
    [Tree, pred] = graphminspantree(sparse(distanceMatrix),StartPoint,'Method','Kruskal');
    % full将稀疏矩阵变为满矩阵，用0填充
    T=full(Tree);
    % 根的前任点还是它自己，从0改为1
    pred(StartPoint)=StartPoint;
    
    % graphconncomp函数，寻找强弱联通子图，程序中的图是无向的
    % S是有多少强连通子图1*1 ， C是这些节点分别属于哪个强连通子图 num(S)*438
    [S, C] = graphconncomp(Tree,'Directed',false);

    % strcat用于连接字符串
    % num2str用于把数字转化成字符串
    details1=strcat('Strongly connected components: ',num2str(S));
    fprintf('%s\n',details1)
    % 统计一下 值 频数 百分比
    tabulate(C)
    
    % path是一个1*438的cell
    path=cell(1,size(distanceMatrix,1));
    % dist是一个1*438的double
    dist=zeros(1,size(distanceMatrix,1));

    for i=1:size(distanceMatrix,1)
        % 遍历每一个节点
        currentI=i;
        path{i}=i;
        % 当当前节点的前任点不是根节点时
        while pred(currentI)~=StartPoint
            % 当前节点的path上加上其前任节点
            path{i}=[path{i},pred(currentI)];
            % 当前节点走到前任节点的位置上
            currentI=pred(currentI);
        end
        % 直到下一步走到根节点
        % 在该path后加上根节点
        path{i}=[path{i},StartPoint];
        
        % dist保存的是从根节点到不同节点的距离，1*438
        dist(i)=0;
        % 对每一段线段
        for j=1:length(path{i})-1
            % index分别表示这个路径中相邻的两个节点
            index1=path{i}(j);
            index2=path{i}(j+1);
            % 不是很懂为什么这里要正反都加一遍，这个T本来就只有下三角部分
            dist(i)=dist(i)+T(index1,index2)+ T(index2,index1);
        end       
    end

    % 装载到结构体属性中
    hopland.dist=dist;
    
    % 如果需要比较
    coef=0;
    if ifdoComparison
        % 算出不同细胞所属阶段 和 根到不同细胞距离 之间的相关系数
        coef = comparison(hopland);
    end 
     
    %% plot tree
    % 打印图    
    if display
        %triplot(dt);
        %hold on
        % 调用了一个自己写的画图程序
        plotMappingResult2DGray(hopland,1);
%         plotMappingResult(hopland,1,0);
        hold on
        %plot
        ms = getoptions(options, 'ms', 5);
        lw = getoptions(options, 'lw', 2);

        %
        for i=1:size(T,1)
            for j=1:size(T,2)
                if T(i,j)~=0
                    path=vertex_orgi(:,[i,j]);
                    h = plot(path(1,:), path(2,:), 'r');
                    set(h,'LineWidth', lw);
                end
            end
        end
    end


end


