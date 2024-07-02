function objPlot(v, f, name)
%STLPLOT is an easy way to plot an STL object
%V is the Nx3 array of vertices
%F is the Mx3 array of faces
%NAME is the name of the object, that will be displayed as a title

figure;
[X, Y, Z] = Convert_3Dmodel(v);
v = [X, Y, Z];
object.vertices = v;
size(object.vertices)
object.faces = f;
size(object.faces)

patch(object,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
hold on

V1 = [-10,-12,0; 15,-12,0; 15,12,0; -10,12,0];
f = [1 2 3 4];
patch('Faces',f,'Vertices',V1,'FaceColor','green','FaceAlpha',.25);
grid off
% hold on
set(gca,'xtick',[],'ytick',[],'ztick',[],'xcolor','w','ycolor','w','zcolor','w')
% V1 = [-1200,-700,300;4400,-700,300;4400,1500,300;-1200,1500,300];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V1,'FaceColor','red','FaceAlpha',.25);
% hold on
% 
% V2 = [-1200,-700,800;4400,-700,800;4400,1500,800;-1200,1500,800];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V2,'FaceColor','red','FaceAlpha',.25);
% hold on
% 
% V3 = [-1200,-700,1300;4400,-700,1300;4400,1500,1300;-1200,1500,1300];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V3,'FaceColor','red','FaceAlpha',.25);
% hold on
% 
% V4 = [-1200,-700,1800;4400,-700,1800;4400,1500,1800;-1200,1500,1800];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V4,'FaceColor','red','FaceAlpha',.25);
% hold on
% 
% V5 = [-1200,-700,2300;4400,-700,2300;4400,1500,2300;-1200,1500,2300];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V5,'FaceColor','red','FaceAlpha',.25);
% hold on


% V1 = [50,0,-10;50,50,-10;50,0,25;50,50,25];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V1,'FaceColor','red','FaceAlpha',.25);
% hold on
% % 
% V2 = [40,0,-10;40,50,-10;40,0,25;40,50,25];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V2,'FaceColor','red','FaceAlpha',.25);
% hold on
% % 
% V3 = [29,0,-10;29,50,-10;29,0,25;29,50,25];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V3,'FaceColor','red','FaceAlpha',.25);
% hold on
% % 
% V4 = [4,0,-10;4,50,-10;4,0,25;4,50,25];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V4,'FaceColor','red','FaceAlpha',.25);
% hold on
% 
% V = [25,0,-10;25,50,-10;25,0,25;25,50,25];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V,'FaceColor','red','FaceAlpha',.25);

xlabel('x')
ylabel('y')
zlabel('z')

% hold on
% V = [36,0,-10;36,50,-10;36,0,25;36,50,25];
% % V5 = [-200,-60,100;-100,-60,100;-200,250,100;-100,250,100];
% f = [1 2 4 3];
% patch('Faces',f,'Vertices',V,'FaceColor','red','FaceAlpha',.25);

% V1 = [-20,0,55;220,0,55;220,140,55;-20,140,55];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V1,'FaceColor','green','FaceAlpha',.2);
% hold on
% V2 = [-20,0,20;220,0,20;220,140,20;-20,140,20];
% f = [1 2 3 4];
% patch('Faces',f,'Vertices',V2,'FaceColor','red','FaceAlpha',.25);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);
grid on;
title(name);