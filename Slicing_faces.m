%% 2023.01.01 by WR
% 该脚本用来实现模型多切片展示
% Input：每个切片之间间隔 △z
% Output：剖面构成的三维图

% [Vertices1, triangles1] = loadOFF_seq('Model4.off');
% V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
% F = triangles1;
% Min_z = min(V(:, 3));
% Max_z = max(V(:, 3));
% Z_set = (Min_z+0.1):1.0:Max_z;
% polygon_set1 = cell(1, length(Z_set));
% polygon_set2 = cell(1, length(Z_set));
% tic;
% parfor i = 1:length(Z_set)
%     [inter_V_1_3D, inter_V_2_3D] = Plot_polygons(Z_set(i));
%     polygon_set1{i} = inter_V_1_3D(:,:);
%     polygon_set2{i} = inter_V_2_3D(:,:);
% end
% toc;
% disp(['计算时间',num2str(toc)]);
% 
% figure
% % hold on
% for ii = 1:length(Z_set)
%     intersect_V_1 = [];
%     inter_V_1 = [];
%     inter_V_2 = [];
%     inter_V_1 = polygon_set1{ii};
%     inter_V_2 = polygon_set2{ii};
%     for iii = 1:length(inter_V_2)
%         intersect_V_1 = [inter_V_1(iii,:);inter_V_2(iii,:)];
%         plot3(intersect_V_1(:,1),intersect_V_1(:,2), intersect_V_1(:,3), 'k', 'LineWidth', 3);
%         hold on
%     end
%     hold on
% end

% x=[76, 127, 190, 380, 759, 1896, 3791];
% a=[0.55063, 0.68203, 0.90375, 1.7162, 3.1738, 7.5006, 14.8861];
% plot(x,a,'--pb',...
%     'LineWidth',4,...
%     'MarkerSize',15,...
%     'MarkerEdgeColor','b')
% 
% % axis([0,8,0,700])  %确定x轴与y轴框图大小
% % set(gca,'XTick',[0:1:8]) %x轴范围1-6，间隔1
% % set(gca,'YTick',[0:100:700]) %y轴范围0-700，间隔100
% % legend('算法1','算法2','算法3');   %右上角标注
% xlabel('The number of cutting faces')  %x轴坐标描述
% ylabel('Time(s)') %y轴坐标描述
% 
% 切片与准确率关系曲线
x=[38, 76, 127, 190, 380, 759, 1896, 3791, 7581, 18951, 37901, 75801, 189502, 379330];
% a=[0.2277, 0.5223, 0.8193, 0.8193, 0.9134, 0.9554, 0.9678, 0.9752, 0.9851, 0.9901, 0.9926, 0.9926, 0.9950, 1];
a=[0.18174, 0.26719, 0.38558, 0.50147, 0.94038, 1.76140, 4.1690, 8.1431, 14.7181, 39.1846, 79.9833, 159.0183, 387.8169, 761.4147];
b=[0.60126, 1.96831, 3.125, 5.160, 12.230, 25.8686, 61.881, 126.2519, 258.4269, 593.8854, 1012.5017, 2099.2167, 4129.788, 8875.20];
c=[0.7349, 2.5563, 4.0016, 5.738, 13.68203, 27.5006, 65.27135, 135.8861, 271.90375, 623.7162, 1236.2519, 2569.3467, 4979.2478, 9563.3753];
% plot(x,a,'-*',...
%     'LineWidth',5,...
%     'MarkerSize',18,...
%     'Color','m',...
%     'MarkerEdgeColor','c')

plot(x,a,'-o',...
    'LineWidth',5,...
    'MarkerSize',18,...
    'Color','r',...
    'MarkerEdgeColor','r')
hold on

plot(x,b,'-x',...
    'LineWidth',5,...
    'MarkerSize',18,...
    'Color','b',...
    'MarkerEdgeColor','b')
hold on

plot(x,c,'-^',...
    'LineWidth',5,...
    'MarkerSize',18,...
    'Color','m',...
    'MarkerEdgeColor','m')

% axis([0,8,0,700])  %确定x轴与y轴框图大小
% set(gca,'xtick',0:100:2500)  
% set(gca,'XTick',[0:10:379003]) %x轴范围1-6，间隔1
% set(gca,'YTick',[0:100:700]) %y轴范围0-700，间隔100
legend('Ours','CuraEngine','FDM');   %右上角标注
set(gca,'FontSize',43);
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlabel('The number of cutting faces')  %x轴坐标描述
% xlabel('The number of slicing faces(log10)')  %x轴坐标描述
% ylabel('Defect detection rate') %y轴坐标描述
ylabel('Time (s)') %y轴坐标描述