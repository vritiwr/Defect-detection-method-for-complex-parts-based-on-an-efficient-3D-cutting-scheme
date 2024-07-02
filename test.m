% tic;
% X=randperm(900);
% X=gpuArray(X);
% % disp(['Before Sort:X=',num2str(X)]);
% disp('--------------------');
% y=for_GPU_vs_CPU(X);
% % disp(['Bubble Sort:x=',num2str(y)]);
% toc;

% tic;
% t = 0.05:0.5:2*pi;
% x1 = cos(t);
% y1 = sin(t);
% x2 = 0.5*cos(t);
% y2 = 0.5*sin(t);
% polyin = polyshape({x1,x2},{y1,y2});
% T = triangulation(polyin);
% toc;

[V,F,VN] = loadOBJ('curve_small1_Defects.obj');
ptCloud = pointCloud(V);

pcwrite(ptCloud,'curve_small1_Defects.ply')