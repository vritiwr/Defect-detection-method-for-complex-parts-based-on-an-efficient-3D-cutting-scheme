%% 2023.03.19 by WR 
%  该脚本主要实现缺陷点处的高亮线段绘制

[Vertices1, triangles1] = loadOFF_seq('Model4.off');
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F = triangles1;
Min_z = min(V(:, 3));
Max_z = max(V(:, 3));
Z_set = (Min_z+0.1):1.0:Max_z;

% [inter_V_1, inter_V_2, Crossed_F] = slicing('Model4.off', 0.4);
% inter_V_1_3D = zeros(length(inter_V_1), 3);
% inter_V_1_3D(:, 1:2) = inter_V_1(:,:);
% inter_V_1_3D(:, 3) = 1;
%        
% inter_V_2_3D = zeros(length(inter_V_2), 3);
% inter_V_2_3D(:, 1:2) = inter_V_2(:,:);
% inter_V_2_3D(:, 3) = 1;
% 
% [inter_V_1_Normal, inter_V_2_Normal, Crossed_F_Normal] = slicing('block.off', 0.4);

for k = 1:length(Z_set)
    [inter_V_1, inter_V_2, Crossed_F] = slicing('Model4.off', Z_set(k));
    [inter_V_1_Normal, inter_V_2_Normal, Crossed_F_Normal] = slicing('block.off', Z_set(k));
    Defect_points_1 = [];
    Defect_points_2 = [];
    Z_z = [Z_set(k); Z_set(k)];
    for i = 1:length(inter_V_1(:,1))
         num1 = 0;
         num2 = 0;
         for j = 1:length(inter_V_1_Normal(:,1))
%         Defect_points = [Defect_points; inter_V_2(j,:)];
            if inter_V_1(i,1) == inter_V_1_Normal(j,1) && inter_V_1(i,2) == inter_V_1_Normal(j,2)
                num1 = num1 + 1;
            end
            if inter_V_2(i,1) == inter_V_2_Normal(j,1) && inter_V_2(i,2) == inter_V_2_Normal(j,2)
                num2 = num2 + 1;
            end
         end
        if num1 == 0
            Defect_points_1 = [Defect_points_1; inter_V_1(i,:)];
        end
        if num2 == 0
            Defect_points_2 = [Defect_points_2; inter_V_2(i,:)];
        end
    end
    if ~isempty(Defect_points_1)
        for ii = 1:length(inter_V_1)
            intersect_V_1 = [inter_V_1(ii,:);inter_V_2(ii,:)];
            plot3(intersect_V_1(:,1),intersect_V_1(:,2), Z_z(:), 'k', 'LineWidth', 3);
            hold on
        end
        hold on
        for jj = 1:length(Defect_points_1)
            intersect_V_defect = [Defect_points_1(jj,:);Defect_points_2(jj,:)];
            plot3(intersect_V_defect(:,1),intersect_V_defect(:,2), Z_z(:), 'r', 'LineWidth', 3);
            hold on
        end
        hold on
    end
end

% 得到缺陷零件的剖面点
Defect_points_1 = [];
Defect_points_2 = [];
for i = 1:length(inter_V_1(:,1))
    num1 = 0;
    num2 = 0;
    for j = 1:length(inter_V_1_Normal(:,1))
%         Defect_points = [Defect_points; inter_V_2(j,:)];
        if inter_V_1(i,1) == inter_V_1_Normal(j,1) && inter_V_1(i,2) == inter_V_1_Normal(j,2)
            num1 = num1 + 1;
        end
        if inter_V_2(i,1) == inter_V_2_Normal(j,1) && inter_V_2(i,2) == inter_V_2_Normal(j,2)
            num2 = num2 + 1;
        end
    end
    if num1 == 0
        Defect_points_1 = [Defect_points_1; inter_V_1(i,:)];
    end
    if num2 == 0
        Defect_points_2 = [Defect_points_2; inter_V_2(i,:)];
    end
end

for ii = 1:length(inter_V_1)
    intersect_V_1 = [inter_V_1(ii,:);inter_V_2(ii,:)];
    plot(intersect_V_1(:,1),intersect_V_1(:,2), 'k', 'LineWidth', 3);
    hold on
end
hold on
for ii = 1:length(Defect_points_1)
    intersect_V_defect = [Defect_points_1(ii,:);Defect_points_2(ii,:)];
    plot(intersect_V_defect(:,1),intersect_V_defect(:,2), 'r', 'LineWidth', 3);
    hold on
end
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Vertices1, triangles1] = loadOFF_seq('Model4.off');
% V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
% F = triangles1;
[V,F] = loadOBJ('4650.obj');
V(:,3) = V(:,3) - 3.9;
Min_z = min(V(:, 3));
Max_z = max(V(:, 3));
Z_set = (Min_z):0.001:Max_z;

figure
for i = 1:length(Z_set)
    [inter_V_1, inter_V_2] = slicing('4650.obj', Z_set(i));
    inter_V_1_3D = zeros(length(inter_V_1), 3);
    inter_V_1_3D(:, 1:2) = inter_V_1(:,:);
    inter_V_1_3D(:, 3) = Z_set(i);

    inter_V_2_3D = zeros(length(inter_V_2), 3);
    inter_V_2_3D(:, 1:2) = inter_V_2(:,:);
    inter_V_2_3D(:, 3) = Z_set(i);
%% show the polygon 
    intersect_V_1 = [];
%   hold on
    for iii = 1:length(inter_V_2)
        intersect_V_1 = [inter_V_1_3D(iii,:);inter_V_2_3D(iii,:)]; 
        plot3(intersect_V_1(:,1),intersect_V_1(:,2),intersect_V_1(:,3), 'Color', '#4DBEEE', 'LineWidth', 3);
        hold on
    end
    hold on

    intersect_V_2 = [];
    [inter_V_1_Normal, inter_V_2_Normal] = slicing('4650.obj', Z_set(i));
    inter_V_1_Normal_3D = zeros(length(inter_V_1_Normal), 3);
    inter_V_1_Normal_3D(:, 1:2) = inter_V_1_Normal(:,:);
    inter_V_1_Normal_3D(:, 3) = Z_set(i);

    inter_V_2_Normal_3D = zeros(length(inter_V_2_Normal), 3);
    inter_V_2_Normal_3D(:, 1:2) = inter_V_2_Normal(:,:);
    inter_V_2_Normal_3D(:, 3) = Z_set(i);

    for ii = 1:length(inter_V_2_Normal)
        intersect_V_2 = [inter_V_1_Normal_3D(ii,:);inter_V_2_Normal_3D(ii,:)];
        plot3(intersect_V_2(:,1),intersect_V_2(:,2),intersect_V_2(:,3), 'Color', '#4DBEEE', 'LineWidth',3);
        hold on
    end
    set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
    set(gca,'linewidth',3.5)
end