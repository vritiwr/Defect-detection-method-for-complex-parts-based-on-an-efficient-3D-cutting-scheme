%% 2023.03.18 by WR 
%  该脚本用来实现计算缺陷点的深度

[Vertices1, triangles1] = loadOFF_seq('Model4.off');
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F = triangles1;

[Vertices2, triangles2] = loadOFF_seq('/Users/wangrui/Desktop/点云项目/off文件读取与导出/block.off');

V_block = reshape(Vertices2, size(Vertices2, 2), size(Vertices2, 3));

defect = [];
Indx = [];
for i = 1:length(V)
    num = 0;
    for j = 1:length(V_block)
        if V(i, 1) == V_block(j, 1) && V(i, 2) == V_block(j, 2) && V(i, 3) == V_block(j, 3)
            num = num + 1;
        end
    end
    if num == 0
        defect = [defect; V(i, :)];
        Indx = [Indx; i];
    end
end

% 找出所有正确的缺陷点坐标
defect_total = [];
for a = 1:length(defect(:,1))
    if defect(a,1) < 0 && defect(a,2) > 0 && defect(a,3) > -2.2
        defect_total = [defect_total; defect(a,:)];
    end
end

% 计算缺陷点距离标准无缺陷模型的最近点
depth = zeros(length(defect_total(:,1)), 1);
for ii = 1:length(defect_total(:,1))
    distance = zeros(length(V_block), 1);
    for jj = 1:length(V_block)
        x1 = defect_total(ii, 1);
        y1 = defect_total(ii, 2);
        z1 = defect_total(ii, 3);
        x2 = V_block(jj, 1);
        y2 = V_block(jj, 2);
        z2 = V_block(jj, 3);
        distance(jj) = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    end
    depth(ii) = min(distance);
end
Max_depth = max(depth);