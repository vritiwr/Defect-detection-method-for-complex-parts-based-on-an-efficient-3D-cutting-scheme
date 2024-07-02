function [inter_V_1_3D, inter_V_2_3D, Crossed_F] = Plot_polygons(z)
[inter_V_1, inter_V_2, Crossed_F] = slicing('Model4.off', z);
inter_V_1_3D = zeros(length(inter_V_1), 3);
inter_V_1_3D(:, 1:2) = inter_V_1(:,:);
inter_V_1_3D(:, 3) = z;

inter_V_2_3D = zeros(length(inter_V_2), 3);
inter_V_2_3D(:, 1:2) = inter_V_2(:,:);
inter_V_2_3D(:, 3) = z;
%% show the polygon


intersect_V_2 = [];
[inter_V_1_Normal, inter_V_2_Normal] = slicing('block.off', z);
for ii = 1:length(inter_V_2_Normal)
    intersect_V_2 = [inter_V_1_Normal(ii,:);inter_V_2_Normal(ii,:)];
    plot(intersect_V_2(:,1),intersect_V_2(:,2), 'k', 'LineWidth', 3);
    hold on
end
hold on
end