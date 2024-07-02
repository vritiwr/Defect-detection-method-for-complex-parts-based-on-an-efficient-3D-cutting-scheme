%% 2023.03.18 by WR
%  该脚本用来实现零件缺陷面积的计算

% [inter_V_1, inter_V_2, Crossed_F] = slicing('Model4.off', 0.1);
% inter_V_1_3D = zeros(length(inter_V_1), 3);
% inter_V_1_3D(:, 1:2) = inter_V_1(:,:);
% inter_V_1_3D(:, 3) = 1;
% 
% inter_V_2_3D = zeros(length(inter_V_2), 3);
% inter_V_2_3D(:, 1:2) = inter_V_2(:,:);
% inter_V_2_3D(:, 3) = 1;
% 
% [inter_V_1_Normal, inter_V_2_Normal, Crossed_F_Normal] = slicing('block.off', 0.2);
% 
% % 得到缺陷零件的缺陷三角形
% Defect_faces = [];
% for i = 1:length(Crossed_F(:,1))
%     num = 0;
%     for j = 1:length(Crossed_F_Normal(:,1))
%         if Crossed_F(i,1) == Crossed_F_Normal(j,1) && Crossed_F(i,2) == Crossed_F_Normal(j,2) && Crossed_F(i,3) == Crossed_F_Normal(j,3)
%             num = num + 1;
%         end
%     end
%     if num == 0
%         Defect_faces = [Defect_faces; Crossed_F(i,:)];
%     end
% end

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

% 找出所有的缺陷三角面片
defect_faces = [];
for ii = 1:length(F(:,1))
    count = 0;
    for jj = 1:length(defect_total(:,1))
        if V(F(ii, 1), 1) == defect_total(jj, 1) && V(F(ii, 1), 2) == defect_total(jj, 2) && V(F(ii, 1), 3) == defect_total(jj, 3)
            count = count + 1;
        end
        if V(F(ii, 2), 1) == defect_total(jj, 1) && V(F(ii, 2), 2) == defect_total(jj, 2) && V(F(ii, 2), 3) == defect_total(jj, 3)
            count = count + 1;
        end
        if V(F(ii, 3), 1) == defect_total(jj, 1) && V(F(ii, 3), 2) == defect_total(jj, 2) && V(F(ii, 3), 3) == defect_total(jj, 3)
            count = count + 1;
        end
    end
    if count >= 2
        defect_faces = [defect_faces; F(ii, :)];
    end
end

Areas = 0;
for aa = 1:length(defect_faces(:,1))
    l1 = 0;
    l2 = 0;
    l3 = 0;
    v1 = defect_faces(aa, 1);
    v2 = defect_faces(aa, 2);
    v3 = defect_faces(aa, 3);
    x1 = V(v1, 1);
    y1 = V(v1, 2);
    z1 = V(v1, 3);
    x2 = V(v2, 1);
    y2 = V(v2, 2);
    z2 = V(v2, 3);
    x3 = V(v3, 1);
    y3 = V(v3, 2);
    z3 = V(v3, 3);
    l1 = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    l2 = sqrt((x1 - x3)^2 + (y1 - y3)^2 + (z1 - z3)^2);
    l3 = sqrt((x2 - x3)^2 + (y2 - y3)^2 + (z2 - z3)^2);
    p = (l1 + l2 + l3) / 2;
    S = sqrt(p * (p - l1) * (p - l2) * (p - l3));
    Areas = Areas + S;
end
