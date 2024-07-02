%% 2023.01.06 by WR
%  该脚本主要实现不同切片情况下的缺陷检出率

[Vertices1, triangles1] = loadOFF_seq('Model4.off');  % defect parts
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));

[Vertices2, triangles2] = loadOFF_seq('block.off');
V_block = reshape(Vertices2, size(Vertices2, 2), size(Vertices2, 3));

Defect_Vertices = [];
V_Indx = [];
for i = 1:length(V)
    num = 0;
    for j = 1:length(V_block)
        if V(i, 1) == V_block(j, 1) && V(i, 2) == V_block(j, 2) && V(i, 3) == V_block(j, 3)
            num = num + 1;
        end
    end
    if num == 0
        Defect_Vertices = [Defect_Vertices; V(i, :)];
        V_Indx = [V_Indx; i];
    end
end

Defect_Faces = [];
Defect_F_Index = [];
for a = 1:length(triangles1)
    Num_V = 0;
    for b = 1:length(V_Indx)
        if triangles1(a,1) == V_Indx(b)
            Num_V = Num_V + 1;
        end
        if triangles1(a,2) == V_Indx(b)
            Num_V = Num_V + 1;
        end
        if triangles1(a,3) == V_Indx(b)
            Num_V = Num_V + 1;
        end
    end
    if Num_V >= 2
        Defect_Faces = [Defect_Faces; triangles1(a,:)];
        Defect_F_Index = [Defect_F_Index; a];
    end
end

Min_z = min(V(:, 3));
Max_z = max(V(:, 3));
Z_set = (Min_z+0.1):0.001:Max_z;
polygon_set1 = cell(1, length(Z_set));
polygon_set2 = cell(1, length(Z_set));
Crossed_Face_set = [];
parfor i = 1:length(Z_set)
    [inter_V_1_3D, inter_V_2_3D, Crossed_F] = Plot_polygons(Z_set(i));
    Crossed_Face_set = [Crossed_Face_set; Crossed_F];
    polygon_set1{i} = inter_V_1_3D(:,:);
    polygon_set2{i} = inter_V_2_3D(:,:);
end
Crossed_Face_set = unique(Crossed_Face_set, 'row');

Faces_with_Defect = [];
for m = 1:length(Crossed_Face_set)
    Num_Defect = 0;
    for n = 1:length(Defect_Faces)
        if Crossed_Face_set(m,1) == Defect_Faces(n,1) && Crossed_Face_set(m,2) == Defect_Faces(n,2) && Crossed_Face_set(m,3) == Defect_Faces(n,3)
            Num_Defect = Num_Defect + 1;
        end
    end
    if Num_Defect == 1
        Faces_with_Defect = [Faces_with_Defect; Crossed_Face_set(m,:)];
    end
end

Defect_rate = length(Faces_with_Defect) ./ length(Defect_Faces);
