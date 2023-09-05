function [inter_V_1, inter_V_2, Crossed_F] = slicing(filename, delta_z)
%% 1.导入文件并将模型旋转、平移
[Vertices1, triangles1] = loadOFF_seq(filename);
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F = triangles1;

% F = F(2:289555,:);

% [M] = max_edge(V,F);            %计算三维模型的最长线段长度
V_N = length(V);                     %计算顶点V的长度
N_F = length(F);                     %计算面片F的维度
[X3, Y3, Z3] = Convert_3Dmodel_integrate(V, delta_z);
V_label = [X3, Y3, Z3];          %将平移后的模型定义一个新矩阵

% figure(1)
% plot3(X3, Y3, Z3);
% hold on
% % % % vvv = [-10,50,0;-10,0,0;70,0,0;70,50,0];
% % % % fff = [1 2 3 4];
% % % % patch('Faces',fff,'Vertices',vvv,'FaceColor','red','FaceAlpha',.3);

N = length(V_label(:,1));        %计算新矩阵的长度
index_v = 1:1:N;                 %构建顶点V的一维标签向量
index_v = index_v';              %将标签向量转置
V_label = [V_label,index_v];     %合并数组得到带有标签的顶点列表
V_original = [X3, Y3, Z3];       %将经过旋转、平移操作后的模型坐标重新赋值给V_original
% figure
xlabel('x')
ylabel('y')
zlabel('z')
name1 = 'cutting model';
% objPlot(V_original, F, name1);
% hold on
% tic;
% m1=cputime;
% % % %% 2.筛选Z轴在[-M,M]之间的点和三角形
% % % % 2.1 筛选点
V_sign = V_label(:,3) > 0;

% % % % 2.2 筛选三角形
% % % new_F = F(unique(new_F_index),:);
% t1 = (cputime-m1)

% m2 = cputime;
%% remove non-crossed triangles  
%% use CPU
Crossed_F = F((abs(sum(V_sign(F),2)-1.5) == 0.5),:);
new_F = F((abs(sum(V_sign(F),2)-1.5) == 0.5),:);
num_F_C = length(new_F);

% t2 = (cputime-m2)
% m3 = cputime;
%% convert each triangle into 2 crossed edges
% newFv_flag=V_original(new_F,3)>0;
%% use CPU
re_new_F = reshape(new_F(:,1:3),[numel(new_F(:,1:3)), 1]); 
new_F_v_z = V_original(re_new_F,3);
new_F_v_z = reshape(new_F_v_z,size(new_F(:,1:3)));
newFv_flag = new_F_v_z > 0;
sub_column = zeros(length(newFv_flag(:,1)),1);
new_F = [new_F,sub_column];

% % if v1/v2 on same side, put v2 to the 4-th column
%% use CPU 
tmp_f = newFv_flag(:,1)==newFv_flag(:,2);
new_F(tmp_f,4) = new_F(tmp_f,2);
new_F(tmp_f,2) = new_F(tmp_f,3);

% % if v1/v3 on same side, put v2 to the 4-th column
%% use CPU
tmp_f = newFv_flag(:,1)==newFv_flag(:,3);
new_F(tmp_f,4) = new_F(tmp_f,2);

% % if v2/v3 on same side, put v2 to the 4-th column
%% use CPU
tmp_f = newFv_flag(:,2)==newFv_flag(:,3);
new_F(tmp_f,4) = new_F(tmp_f,1);

% t3 = (cputime-m3)
% m4 = cputime;

%% compute the intersect points
%% use CPU
edgeLen=abs(V_original(new_F(:,1),3)) + abs(V_original(new_F(:,2),3));
inter_V_1 = V_original(new_F(:,2),[1 2]) .* [abs(V_original(new_F(:,1),3))./edgeLen abs(V_original(new_F(:,1),3))./edgeLen] ...
			+ V_original(new_F(:,1),[1 2]) .* [abs(V_original(new_F(:,2),3))./edgeLen abs(V_original(new_F(:,2),3))./edgeLen];

edgeLen=abs(V_original(new_F(:,3),3)) + abs(V_original(new_F(:,4),3));
inter_V_2 = V_original(new_F(:,4),[1 2]) .* [abs(V_original(new_F(:,3),3))./edgeLen abs(V_original(new_F(:,3),3))./edgeLen] ...
			+ V_original(new_F(:,3),[1 2]) .* [abs(V_original(new_F(:,4),3))./edgeLen abs(V_original(new_F(:,4),3))./edgeLen];

% t4 = (cputime-m4)

% toc;
% disp(['计算时间',num2str(toc)]);

end