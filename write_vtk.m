function write_vtk(V,F,n,i)
% write_vtk (V,F,n,i) saves meshes to disk in vtk format. 
% V is the matrix of vertices 
% F is the matrix of faces 
% n is a number that help us to create the file name. It could be {1,2,3} 
    % 1 = left hippocampus
    % 2 = right hippocampus
    % 3 = aligned mesh

% i is the iteration that helps us to maintain a relationship 
% with the source file
 
% We calculate the size of vertices and faces matrix (in this case only
% the number of columns 
[~,n1]=size(V);
[~,n2]=size(F);


% In order to write the face matrix into the file, we have to write in the
% firt column the number of elements that form one face. Like the meshes
% are created with triangles, we only need 3 points. So we have to write
% the number 3. 
tr=3*ones(1,n2); %vector of 3 with the same dimensions of faces. 
F_men= F - 1; %reduce one unit each position         
F_g=[tr; F_men]; %add the vector of 3's


% The name of de file depends on the variable n. 
if n == 1 %it order to save a left hippocampus
    filename= sprintf('mesh%d_left.vtk',i);
end 
if n == 2 %right hippocampus
    filename = sprintf('mesh%d_right.vtk',i);
end
if n == 3 %an aligned mesh
    filename= sprintf('aligned_mesh%d.vtk',i); 
end 
    

% First, we have to open the file. If the file doesn't exist, it will
% create one.
file = fopen(filename, 'w');
if file == -1
    error('Cannot open file for writing.');
    return;
end

% It writes a header with useful information like the format vtk or the
% type of code: ASCII 
fprintf(file,'# vtk DataFile Version 4.0\n');
fprintf(file,'vtk output\n');
fprintf(file,'ASCII\n');
fprintf(file,'DATASET POLYDATA\n');

% To write the vertice matrix
fprintf(file,'POINTS %d float \n',n1);
fprintf(file,'%f %f %f\n', V(:));

%To write the face matrix
fprintf(file,'POLYGONS %d %d\n', n2, 4*n2);
fprintf(file,'%d %d %d %d\n', F_g(:));

%Finally, it closes the file.
fclose(file);
fprintf(' Save it!\n')


end
