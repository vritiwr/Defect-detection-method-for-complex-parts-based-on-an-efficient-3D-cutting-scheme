

[Vertices1, triangles1] = loadOFF_seq('Model4.off');
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F = triangles1;
% [V,F] = loadOBJ('female_0000.obj');
% % % data = struct('vertices',V,'faces',F);


%Save the mesh to a new stl file:
% WRITE_stl('sample2.stl',coordVERTICESnew,coordNORMALSnew)
% 
% %Plot the mesh with the normals:
% % figure;
% % PLOT_3D_stl_patch(coordVERTICESnew,coordNORMALSnew)
% 
%Write 3D mesh to a new ply file
ptCloud = pointCloud(V);
pcwrite(ptCloud,'synthesis.ply')
ptObj = pcread('synthesis.ply');

% Compute FPFH feature
keyInds = 1:1:2324;
% keyInds = 1:1:2;
features = extractFPFHFeatures(ptObj,keyInds);



% extract defect points
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

% USing FCM algorithm to divide the data into 2 categroies
% [Umat,Cmat,it,valJ] = fcm(features,2,2,0.0001,100)
metric = @euclidean;
Max = 1000;
tol = 1e-3;
[prediction,v] = fcm_1(2, features, 2, metric, Max, tol);

% caculate the defected points by FPFH and FCM
defect_points = [];
normal_points = [];
for i = 1:length(V)
    if prediction(i) == 2
        defect_points = [defect_points; V(i, :)];
    else
        normal_points = [normal_points; V(i, :)];
    end
end

defect_total = [];
for a = 1:length(defect(:,1))
    if defect(a,1) < 0 && defect(a,2) > 0 && defect(a,3) > -2.2
        defect_total = [defect_total; defect(a,:)];
    end
end

% map = addcolorplus(300);
figure
scatter3(V(:,1), V(:,2), V(:,3), 'MarkerEdgeColor', [0.26953125 0.53515625 0.578125], 'MarkerFaceColor',[0.26953125 0.53515625 0.578125])
% scatter3(V(:,1), V(:,2), V(:,3), 'MarkerEdgeColor', [0.5171875 0.68359375 0.60546875], 'MarkerFaceColor',[0.5171875 0.68359375 0.60546875])
% scatter3(V(:,1), V(:,2), V(:,3), 'MarkerEdgeColor', [0 .75 .75], 'MarkerFaceColor',[0 .75 .75])
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
hold on
scatter3(defect_total(:,1), defect_total(:,2), defect_total(:,3), 'MarkerEdgeColor', [0.9921875 0.26172 0.39453125],'MarkerFaceColor',[0.9921875 0.26172 0.39453125])
hold on
% scatter3(defect_points(:,1), defect_points(:,2), defect_points(:,3), 'MarkerEdgeColor', [0.97265625 .80078125 .67578125], 'MarkerFaceColor',[0.97265625 .80078125 .67578125])
