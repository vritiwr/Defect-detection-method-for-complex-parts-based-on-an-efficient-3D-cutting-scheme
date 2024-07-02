%% 2019.12.09 by 王睿
%% 主要功能为读出stl文件中的数据

[f,n,v]=stlread('01.bottle.stl');
F = f';
V = v';
patch(f,n,v)