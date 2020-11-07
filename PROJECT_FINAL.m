%% TASK 1

%We create a cell 'M' with all the meshes that contain both hippocampus

M={'both_1.vtk','both_2.vtk','both_3.vtk', 'both_4.vtk', 'both_5.vtk',  ...
     'both_6.vtk', 'both_7.vtk','both_8.vtk', 'both_9.vtk', 'both_10.vtk', ...
     'both_11.vtk', 'both_12.vtk', 'both_13.vtk', 'both_14.vtk', 'both_15.vtk',...
     'both_16.vtk', 'both_17.vtk', 'both_18.vtk', 'both_19.vtk', 'both_20.vtk'}; 

%In order to see all the meshes with both hippocampus that we have we 
%create two empty cells
DiccV = cell(20,1,1);
DiccF = cell(20,1,1);

figure;
for i = 1:length(M) %it goes through all the cell 'M' 
     %Loads the corresponding mesh
    [DiccV{i}, DiccF{i}] =read_vtk(M{i}); %read data from VTK file
    
    
    %To plot all the meshes
    sgtitle('All the meshes')
    subplot(5,5,i)
    plt_mesh(DiccV{i},DiccF{i}); 
    hold on;
    
    
end

%% TASK 1
%We want to separate both hippocampus, so we iterate through all the 
%cell 'M' in order to load the meshes and work with them.  
for i = 1:length(M)
    mesh = M{i}; %we save one of the meshes in a variable 
    [V, F] =read_vtk(mesh); %read data from VTK file

% Vertices, array size n.vertices x 3. Each row represents a vertice, each
% column is the x,y or z coordinate in 3D space.
whos V

% Faces, array size n.vertices x3 specifying the connectivity of the mesh.
% Each row represents a triangle, and the three columns specify the three
% vertices that connect to it
whos F

% Visualize the mesh
options.edge_color=0;

%Plot the mesh 
figure;
plt_mesh(V,F,options)
str = sprintf('Both hippocampus of mesh %d', i);
title(str)


%To find the connected components in the mesh
[F1,V1,F2,V2] = find_cc(V,F);


% V1 and F1 correspond to the left hippocampus 
figure;
plt_mesh(V1,F1,options) %plot the left part
str_left = sprintf('Left hippocampus of mesh %d', i);
title(str_left)


% V2 and F2 correspond to the right hippocampus
figure;
plt_mesh(V2,F2,options) %plot the right part
str_right = sprintf('Right hippocampus of mesh %d', i);
title(str_right)


%TASK 1 - OPTIONAL - To save vtk meshes to disk
% If we want to save a mesh from the faces that we obtained in find_cc, 
%we have transpose the Face matrix. 

write_vtk(V1,F1',1,i) %save left hippocampus mesh into vtk file

write_vtk(V2,F2',2,i) %save the right one

end 

    
%% TASK 2 - LEFT PART

%Load two different meshes that only contain one hippocampus in this case
%the left one. 
mesh1 = 'mesh17_left.vtk';
mesh2 = 'mesh18_left.vtk';

%Read both hippocampus and save it into 2 matrices: Vetices and Faces
[VL1, FL1] =read_vtk(mesh1); 
[VL2, FL2] =read_vtk(mesh2);

whos VL1; whos VL2;

whos FL1; whos FL2;

% Visualize the mesh
options.edge_color=0;

figure;
subplot (121)
plt_mesh(VL1,FL1,options)
title('Left hippocampus - REFERENCE mesh nº 17')
subplot (122)
plt_mesh(VL2,FL2,options)
title('Left hippocampus - mesh nº 18')

[~,VL3] = procrustes(VL1',VL2'); %The first input is the reference, the second is the mesh to align.
%For the output we only need the vertex matrix, because the faces will
%remain the same, the only thing that will change is the location of the
%vertex.
figure;
plt_mesh(VL3,FL1,options); %Alineation V2 with V1 -> V3 are the new aligned vertices . The faces do not change.
title('Aligned left hippocampus') 

%To see if it the alineation is done correctly
figure;
plt_mesh(VL1,FL1,options) %Original one
hold on;
plt_mesh(VL3,FL1,options) %After alineation

%Save it in a vtk file
%write_vtk(VL3',FL1,1,221);

%% TASK 2 - RIGHT PART

%Load two right hippocampus meshes
mesh3 = 'mesh17_right.vtk';
mesh4 = 'mesh18_right.vtk';

%Read them 
[VR1, FR1] =read_vtk(mesh3);
[VR2, FR2] =read_vtk(mesh4);


whos VR1; whos VR2;

whos FR1; whos FR2;

% Visualize the mesh 
figure;
subplot (121)
plt_mesh(VR1,FR1,options)
title('Right hippocampus - REFERENCE mesh nº 17')
subplot (122)
plt_mesh(VR2,FR2,options)
title('Right hippocampus - mesh nº 18')

[~,VR3] = procrustes(VR1',VR2'); %The first input is the reference 
figure;
plt_mesh(VR3,FR1,options); %Alineation V2 with V1 -> V3 are the new aligned vertices . The faces do not change.
title('Aligned right hippocampus') 

%To see if it the alineation is done correctly
figure;
plt_mesh(VR2,FR2,options) %original one
hold on;
plt_mesh(VR3,FR1,options) %after alineation

%Save it in a vtk file
%write_vtk(VR3',FR1,2,22);

%% TASK 2 - OPTIONAL
% LEFT PART
[hippo,VL3_O] = newProcrustes(VL1',VL2'); %The first input is the reference 
figure;
plt_mesh(VL3_O,FL1,options); %Alineation V2 with V1 -> V3 are the new aligned vertices . The faces do not change.
title('Aligned left hippocampus - OPTIONAL') 

%To see if it the alineation is done correctly
figure;
plt_mesh(VL2,FL2,options) %original one
plt_mesh(VL3_O,FL1,options) %after alineation

%Save it in a vtk file
%write_vtk(VL3_O',FL1,1,22);


% RIGHT PART
[hippo2,VR3_O] = newProcrustes(VR1',VR2'); %The first input is the reference 
figure;
plt_mesh(VR3_O,FR1,options); %Alineation V2 with V1 -> V3 are the new aligned vertices . The faces do not change.
title('Aligned right hippocampus - OPTIONAL') 

%To see if it the alineation is done correctly
figure;
plt_mesh(VR2,FR2,options) %original one
plt_mesh(VR3_O,FR1,options) %after alineation

%Save it in a vtk file
%write_vtk(VR3_O',FR1,2,22);
%% TASK 3 - LEFT

%1. Choose a mesh as reference.
%Charge all the meshes into a cell.
[Vref, Fref] =read_vtk('mesh3_left.vtk');
M_left={'mesh1_left.vtk','mesh2_left.vtk','mesh3_left.vtk','mesh4_left.vtk',...
    'mesh5_left.vtk','mesh6_left.vtk','mesh7_left.vtk','mesh8_left.vtk',...
    'mesh9_left.vtk','mesh10_left.vtk','mesh11_left.vtk','mesh12_left.vtk',...
    'mesh13_left.vtk','mesh14_left.vtk','mesh15_left.vtk','mesh16_left.vtk',...
    'mesh17_left.vtk','mesh18_left.vtk','mesh19_left.vtk','mesh20_left.vtk'};

Diccoriginals = cell(20,1,1);
Diccaligned = cell(20,1,1); %The aligned meshes will be saved here.
Vsum=0;

%Read all meshes and save it to the originals cell
for i = 1:length(M_left)
    [Diccoriginals{i},~] = read_vtk(M_left{i}); 
end


for k = 1:length(Diccoriginals)
    
    for i = 1: length(Diccoriginals)
        %2. Compute Procrustes analysis of all meshes to reference mesh to align all meshes to the reference.
        %To align all of them and compute the mean shape correctly, 
        [~,Diccaligned{i}] = procrustes(Vref',Diccoriginals{i}');

        %3. Compute mean shape of the aligned meshes.
        VSum = Vsum+Diccaligned{i};
        VMean=Diccaligned{i}/i; %Mean shape computation
    end
   
    %Once the Mean is calculated, the norm of each column of each matrix
    %should be calculated
    x=norm(Vref(1,:)');
    y=norm(Vref(2,:)');
    z=norm(Vref(3,:)');
    
    u=norm(VMean(:,1));
    v=norm(VMean(:,2));
    w=norm(VMean(:,3));  
    
    
    % 4. If the Procrustes distance between the new mean shape and the reference shape is below a set threshold, 
    %stop. Else, use the obtained mean shape as reference and return to step 2 (-> Beggining of the loop).
    
    
    d = sqrt((u-x)^2+(v-y)^2+(w-z)^2); %Distance calulation
        if d>0.3 %Threshold
            %If the distance is above the threshold, the Mean should be
            %recalculated
            Vref=VMean';
        else
            break
        end
       
end

%All the lefts hippocampus superposed 
for j = 1:length(Diccaligned)
    plt_mesh(Diccaligned{j},Fref);
    hold on;
end

%All the lefts hippocampus aligned
figure;
for m = 1:length(Diccaligned)
    sgtitle('Aligned')
    subplot(2,10,m)
    plt_mesh(Diccaligned{m},Fref);
    hold on;
  
    %write_vtk(Diccaligned{m}', Fref,3,m); 
end

%% TASK 3 - RIGHT


%1. Choose a mesh as reference.

[Vref, Fref] =read_vtk('mesh3_right.vtk');
M_right={'mesh1_right.vtk','mesh2_right.vtk','mesh3_right.vtk','mesh4_right.vtk',...
    'mesh5_right.vtk','mesh6_right.vtk','mesh7_right.vtk','mesh8_right.vtk',...
    'mesh9_right.vtk','mesh10_right.vtk','mesh11_right.vtk','mesh12_right.vtk',...
    'mesh13_right.vtk','mesh14_right.vtk','mesh15_right.vtk','mesh16_right.vtk',...
    'mesh17_right.vtk','mesh18_right.vtk','mesh19_right.vtk','mesh20_right.vtk'};

Diccoriginals = cell(20,1,1); 
Diccaligned = cell(20,1,1); 
Vsum=0;

%To read all the meshes
for i = 1:length(M_right)
    [Diccoriginals{i},~] = read_vtk(M_right{i}); 
end


for k = 1:length(Diccoriginals)
    
    for i = 1: length(Diccoriginals)
        %2. Compute Procrustes analysis of all meshes to reference mesh to align all meshes to the reference.
        [~,Diccaligned{i}] = procrustes(Vref',Diccoriginals{i}');

        %3. Compute mean shape of the aligned meshes.
        VSum = Vsum+Diccaligned{i};
        VMean=Diccaligned{i}/i;
    end
   
    
    x=norm(Vref(1,:)');
    y=norm(Vref(2,:)');
    z=norm(Vref(3,:)');
    
    u=norm(VMean(:,1));
    v=norm(VMean(:,2));
    w=norm(VMean(:,3));  
    
    
    % 4. If the Procrustes distance between the new mean shape and the reference shape is below a set threshold, 
    %stop. Else, use the obtained mean shape as reference and return to step 2 (-> Beggining of the loop).
    
    
    d = sqrt((u-x)^2+(v-y)^2+(w-z)^2) %distance calculation
        if d>0.3 %threshold
            Vref=VMean';
        else
            break
        end
       
end

%All the roght hippocampus superposed 
for j = 1:length(Diccaligned)
    plt_mesh(Diccaligned{j},Fref);
    hold on;
end

%All the lefts hippocampus aligned 
figure;
for m = 1:length(Diccaligned)
    sgtitle('Aligned')
    subplot(2,10,m)
    plt_mesh(Diccaligned{m},Fref); 
    hold on;
    
    %write_vtk(Diccaligned{m}', Fref,3,m);
end



