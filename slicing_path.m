%% 2019.07.12 by wr
%% 调用剖切算法对obj文件进行切割  usage：slicing_path('.\slicing\objfiles\')
function slicing_path(fileDir,nbOBJ)
fileDir_V = [fileDir 'V\'];
fileDir_F = [fileDir 'F\'];
files_V = dir([fileDir_V '*.mat']);
files_F = dir([fileDir_F '*.mat']);
% nbOBJs = length(files);
Times = [];
num_V = 0;
num_F = 0;
num_c_tri = 0;
num_c_p = 0;
V_{nbOBJ}=[];
F_{nbOBJ}=[];
M_=zeros(1,nbOBJ);
for i = 1:nbOBJ
%     figure(1)
    V = load([fileDir_V files_V(i).name]);
    V = V.V;
    V_{i}=V;
    F = load([fileDir_F files_F(i).name]);
    F = F.F;
    F_{i}=F;
%     [V,F] = loadOBJ_slicing([fileDir files(i).name]);
%     plot3(V(:,1),V(:,2),V(:,3))
%     hold on
%     vvv = [-20,0,0;140,0,0;140,140,0;-20,140,0];
%     fff = [1 2 3 4];
%     patch('Faces',fff,'Vertices',vvv,'FaceColor','red','FaceAlpha',.3);
    [t1,t2,t3,t4,inter_V_1,inter_V_2,V_N,N_F,num_F_C,points_num] = slicing0731(V,F);
% % %     M_(i)=M;
    Times(i,:) = [t1,t2,t3,t4];
    num_V = num_V + V_N;
    num_F = num_F + N_F;
    num_c_tri = num_c_tri + num_F_C;
    num_c_p = num_c_p + points_num;
    intersect_V = [];
%     figure(2)
%     axis([-10,220,-10,220]);
%     xlim([-10,220]);%对X轴设定显示范围 
    hold on
    for iii = 1:length(inter_V_2)
        intersect_V = [inter_V_1(iii,:);inter_V_2(iii,:)];
        plot(intersect_V(:,1),intersect_V(:,2));
        hold on
    end
    hold on
%     if ~isempty(inter_V_1)
%         tri = delaunay(inter_V_1(:,1),inter_V_1(:,2));
%         hold on
%         triplot(tri,inter_V_1(:,1),inter_V_1(:,2));
%     end
end
T1 = sum(Times(:,1))
T2 = sum(Times(:,2))
T3 = sum(Times(:,3))
T4 = sum(Times(:,4))
% T5 = sum(Times(:,5))
% T5 = sum(Times(:,5))
TT = T1 + T2 + T3 + T4
num_V
num_F
num_c_tri
num_c_p
end