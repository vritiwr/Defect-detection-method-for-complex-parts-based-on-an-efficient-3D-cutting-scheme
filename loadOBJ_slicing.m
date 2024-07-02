function [V,F] = loadOBJ_slicing(filename)
  % Reads a .obj mesh file and outputs the vertex and face list
  % assumes a 3D triangle mesh and ignores everything but:
  % v x y z and f i j k lines
  % Input:
  %  filename  string of obj file's path
  %
  % Output:
  %  V  number of vertices x 3 array of vertex positions
  %  F  number of faces x 3 array of face indices
  %
  %Usage :  loadOBJ('D:\2DPDF\sub_obj\1\3167}_0.obj')   C:\Users\hp\Desktop\slicing\slicing\3167}_0.obj
  
  V = zeros(0,3);
  F = zeros(0,3);
  vertex_index = 1;
  face_index = 1;
  fid = fopen(filename,'rt');
  line = fgets(fid);
%   a = 0;
  while ischar(line)
    vertex = sscanf(line,'v %f %f %f');
%     a = a + 1;    %calculating the number of vertices in the obj file
%     face = sscanf(line,'f %d %d %d');
%     face_long = sscanf(line,'f %d %d %d');
    face_long1 = sscanf(line,'f %d/%d/%d  %d/%d/%d  %d/%d/%d');
    face_long2 = sscanf(line,'f %d//%d  %d//%d  %d//%d');
    
    % see if line is vertex command if so add to vertices
    if(size(vertex)>0)
      V(vertex_index,:) = vertex;
      vertex_index = vertex_index+1;
%     % see if line is simple face command if so add to faces
%     elseif(size(face)>0)
%       F(face_index,:) = face;
%       face_index = face_index+1;
    % see if line is a long face command if so add to faces
    elseif((length(face_long1)>0) || (length(face_long2)>0))
      % remove normal and texture indices
%       face_long = face_long(1:2:end);
%       F(face_index,:) = face_long(1:end);
      if ~isempty(face_long1) && face_index < 81
          F1(face_index,:) = face_long1(3:3:end);
      end
      if ~isempty(face_long2) && face_index > 80
          F2(face_index,:) = face_long2(2:2:end);
      end
      face_index = face_index+1;
    else
      donothing=1;
    end

    line = fgets(fid);
  end
  F = [F1;F2];
  F = [F(1:80,:);F(161:240,:)];
  save('Faces.mat','F');
  fclose(fid);

end