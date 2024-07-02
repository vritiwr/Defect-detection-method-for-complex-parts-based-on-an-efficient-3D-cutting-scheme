%% ICP case
[Vertices1, triangles1] = loadOFF_seq('Model4.off');
V1 = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F1 = triangles1;

[Vertices2, triangles2] = loadOFF_seq('block.off');
V2 = reshape(Vertices2, size(Vertices2, 2), size(Vertices2, 3));
F2 = triangles2;

[R, T] = ICP(V1, V2);
newdata = R*V2' + T;
V = newdata';

figure;
plot3(V1(:,1), V1(:,2), V1(:,3), '.', 'color', '#FEB2B4', 'MarkerSize', 20)
hold on
plot3(V(:,1), V(:,2), V(:,3), '.', 'color', '#9394E7', 'MarkerSize', 20)
hold off
set(gca, 'xtick', [], 'ytick', [], 'ztick', [], 'xcolor', 'w', 'ycolor', 'w', 'zcolor', 'w')

