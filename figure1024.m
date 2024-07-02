% [v, f] = Copy_of_loadOBJ_slicing('building.obj');
% name = 'building.obj';
% objPlot(v, f, name);

% [V, F] = Copy_of_loadOBJ_slicing('Farm_house.obj');
% name1 = 'Farm_house.obj';
% objPlot(V, F, name1);


[Vertices1, triangles1] = loadOFF_seq('Model4.off');
V = reshape(Vertices1, size(Vertices1, 2), size(Vertices1, 3));
F = triangles1;
name1 = 'Farm_house.obj';
objPlot(V, F, name1);

% [V1, F1] = Copy_of_loadOBJ_slicing('windmill.obj');
% name2 = 'windmill.obj';
% objPlot(V1, F1, name2);

% [V2, F2] = Copy_of_loadOBJ_slicing('San Fran Tower.obj');
% name2 = 'San Fran Tower.obj';
% objPlot(V2, F2, name2);

% [V3, F3] = Copy_of_loadOBJ_slicing('part1model.obj');
% name3 = 'part1model.obj';
% objPlot(V3, F3, name3);

% t = 0.05:0.5:2*pi;
% x1 = cos(t);
% y1 = sin(t);
% x2 = 0.5*cos(t);
% y2 = 0.5*sin(t);
% polyin = polyshape({x1,x2},{y1,y2})
% plot(polyin)
% hold on
% scatter(x1,y1,45,'filled','r')
% hold on
% scatter(x2,y2,45,'filled','r')z