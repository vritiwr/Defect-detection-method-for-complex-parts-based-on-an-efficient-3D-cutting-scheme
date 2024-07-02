%% 2023.01.02 by WR
% 该脚本用来计算有缺陷零件与无缺陷零件对比后的缺陷点

[inter_V_1, inter_V_2] = slicing('Model4.off', 0.1);
inter_V_1_3D = zeros(length(inter_V_1), 3);
inter_V_1_3D(:, 1:2) = inter_V_1(:,:);
inter_V_1_3D(:, 3) = 1;

inter_V_2_3D = zeros(length(inter_V_2), 3);
inter_V_2_3D(:, 1:2) = inter_V_2(:,:);
inter_V_2_3D(:, 3) = 1;
%% show the polygon
intersect_V_1 = [];
figure
% hold on
for iii = 1:length(inter_V_2)
    intersect_V_1 = [inter_V_1(iii,:);inter_V_2(iii,:)]; 
    plot(intersect_V_1(:,1),intersect_V_1(:,2), 'r', 'LineWidth', 6);
    hold on
end
hold on

intersect_V_2 = [];
[inter_V_1_Normal, inter_V_2_Normal] = slicing('block.off', 0.1);
for ii = 1:length(inter_V_2_Normal)
    intersect_V_2 = [inter_V_1_Normal(ii,:);inter_V_2_Normal(ii,:)];
    plot(intersect_V_2(:,1),intersect_V_2(:,2), 'Color', '#4DBEEE', 'LineWidth',6);
    hold on
end
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
set(gca,'linewidth',3.5)

% 得到缺陷零件的剖面点
Defect_points = [];
for i = 1:length(inter_V_1(:,1))
    num = 0;
    for j = 1:length(inter_V_2(:,1))
        Defect_points = [Defect_points; inter_V_2(j,:)];
        if inter_V_1(i,1) == inter_V_2(j,1) && inter_V_1(i,2) == inter_V_2(j,2)
            num = num + 1;
        end
    end
    if num == 0
        Defect_points = [Defect_points; inter_V_1(i,:)];
    end
end

Defect_points = unique(Defect_points, 'rows');

% 得到正常零件的剖面点
normal_points = [];
for m = 1:length(inter_V_1_Normal(:,1))
    num_normal = 0;
    for n = 1:length(inter_V_2_Normal(:,1))
        normal_points = [normal_points; inter_V_2_Normal(n,:)];
        if inter_V_1_Normal(m,1) == inter_V_2_Normal(n,1) && inter_V_1_Normal(m,2) == inter_V_2_Normal(n,2)
            num_normal = num_normal + 1;
        end
    end 
    if num_normal == 0
        normal_points = [normal_points; inter_V_1_Normal(m,:)];
    end
end

normal_points = unique(normal_points, 'rows');

% 通过将有无缺陷零件的剖面点对齐的到缺陷部分
Points_with_defects = [];
for a = 1:length(Defect_points(:,1))
    count = 0;
    for b = 1:length(normal_points(:,1))
        if Defect_points(a,1) == normal_points(b,1) && Defect_points(a,2) == normal_points(b,2)
            count = count + 1;
        end
    end
    if count == 0
        Points_with_defects = [Points_with_defects; Defect_points(a,:)];
    end
end

% x = Points_with_defects(:,1);
% y = Points_with_defects(:,2);
% sz = 25;
% % c = linspace(1,10,length(x));
% scatter(x,y,sz,'filled')

points_set1 = [];
for c = 1:length(inter_V_1(:,1))
    for d = 1:length(Points_with_defects(:,1))
        if inter_V_1(c,1) == Points_with_defects(d,1) && inter_V_1(c,2) == Points_with_defects(d,2)
            points_set1 = [points_set1; inter_V_1(c,:)];
        end
    end
end

points_set2 = [];
for c = 1:length(inter_V_2(:,1))
    for d = 1:length(Points_with_defects(:,1))
        if inter_V_2(c,1) == Points_with_defects(d,1) && inter_V_2(c,2) == Points_with_defects(d,2)
            points_set2 = [points_set2; inter_V_2(c,:)];
        end
    end
end

intersect_defect_points = [];
for aa = 1:length(points_set2)
    intersect_defect_points = [points_set1(aa,:);points_set2(aa,:)];
    plot(intersect_defect_points(:,1),intersect_defect_points(:,2), 'b');
    hold on
end