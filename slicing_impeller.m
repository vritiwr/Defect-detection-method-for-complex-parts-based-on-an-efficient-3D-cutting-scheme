%% 2023.03.20 by WR
%  该脚本主要用来实现叶轮零件的剖切结果

% [V,F,VN] = loadOBJ('failure.obj');
% save('failure_vertex.mat', 'V')
% save('failure_face.mat', 'F')
VV = load('failure_vertex.mat');
V = VV.V;
FF = load('failure_face.mat');
F = FF.F;
V(:,3) = V(:,3) - 325.09;
% V(:,3) = V(:,3) - 239;
% ptCloud = pointCloud(V);
% pcwrite(ptCloud,'4.ply')
% ptObj = pcread('4.ply');
% gridStep = 5;
% ptCloudA = pcdownsample(ptCloud,'gridAverage',gridStep);
% V1 = ptCloudA.Location;
% ptCloud = pointCloud(V1);

% pcwrite(ptCloud,'4_downsample.ply')
% ptObj = pcread('4_downsample.ply');

% DT = delaunay(V1(:,1),V1(:,2));
% trisurf(DT,V1(:,1),V1(:,2),V1(:,3))
% TR = triangulation(DT,V1(:,1),V1(:,2),V1(:,3));

Min_z = min(V(:, 3));
Max_z = max(V(:, 3));
Z_ = (Min_z+0.1):0.6:Max_z;
figure
for i = 1:length(Z_)

    intersect_V_2 = [];
    [inter_V_1_Normal, inter_V_2_Normal] = slicing_YeLun(V, F, Z_(i));
    inter_V_1_Normal_3D = zeros(length(inter_V_1_Normal), 3);
    inter_V_1_Normal_3D(:, 1:2) = inter_V_1_Normal(:,:);
    inter_V_1_Normal_3D(:, 3) = Z_(i);

    inter_V_2_Normal_3D = zeros(length(inter_V_2_Normal), 3);
    inter_V_2_Normal_3D(:, 1:2) = inter_V_2_Normal(:,:);
    inter_V_2_Normal_3D(:, 3) = Z_(i);

    for ii = 1:length(inter_V_2_Normal)
        intersect_V_2 = [inter_V_1_Normal_3D(ii,:);inter_V_2_Normal_3D(ii,:)];
        plot3(intersect_V_2(:,1),intersect_V_2(:,2),intersect_V_2(:,3), 'Color', '#4DBEEE', 'LineWidth',3);
        hold on
    end
    set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
    set(gca,'linewidth',3.5)
end