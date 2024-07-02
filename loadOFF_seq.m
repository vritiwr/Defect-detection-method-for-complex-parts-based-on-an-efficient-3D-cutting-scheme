% 2022/10/10
% This function loads a list of '.off' files.
% 
% 
% usage: 
% [Vertices,triangles]=loadOFF_seq(fileDir)

function [Vertices,triangles]=loadOFF_seq(fileDir)

disp('/////////////////////LOADING OFF SEQUENCE/////////////////')
files=dir(fileDir);
% files(1:2)=[];
nbF=length(files);
% [V,triangles] = loadOBJ([fileDir '/' files(1).name]);
% [vertex,face] = read_off(filename)
[V,triangles] = read_off(files(1).name);
V=V';triangles=triangles';
nbV=size(V,1);

Vertices=zeros(nbF,nbV,3);
Vertices(1,:,:)=V;

for ii=2:nbF
% 	animVertices(ii,:,:) = loadOBJ([fileDir '/' files(ii).name]);
    tt = read_off([fileDir '\' files(ii).name]);
    Vertices(ii,:,:)=tt';
end 

% % %% export obj_seq
% % filePath='..\cloth\ORIG\OBJ';
% % exportOBJseq( filePath, animVertices, triangles)