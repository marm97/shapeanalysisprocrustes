function [d, Z] = newProcrustes(X, Y)
% [d, Z] = newProcrustes(X, Y) determines a linear transformation (translation,
%   reflection, orthogonal rotation, and scaling) of the points in the
%   matrix Y to best conform them to the points in the matrix X. 
%   d is the component for which the linear tranformation is done and Z is
%   the vertex output as new aligned mesh. To this vertex, a faces matrix must be added with read_vtk method.
%   X and Y are the two vertex matrices to be aligned, X the reference, and
%   Y the one that will suffer the alignment changes.


%We save the dimentions of the two input matrices, even If the case of nX
%different from nY, and mX different from mY is not contemplated (no padding of zeros)
[nX, mX]   = size(X);   
[~, mY] = size(Y);


%WE CALCULATE THE CENTER (XC,YC):

%Horizontal vectors containing in each column the mean of each column of the original matrix of vectors
columnMeanX = mean(X,1);    
columnMeanY = mean(Y,1);

%repmat replicates the 2 row vectors with the means 3 times, so it is
%possible to subtract the means from the 2 original matrics 3 x numbers of
%columns
XC = X - repmat(columnMeanX, nX, 1);
YC = Y - repmat(columnMeanY, nX, 1);

%Sum of the squares of each column; I get a line vector
squareX = sum(XC.^2,1);
squareY = sum(YC.^2,1);

%Sum of all the values of the row vector, getting a single value
squareX = sum(squareX);
squareY = sum(squareY);

%Square roots of the 2 values
normX = sqrt(squareX);  
normY = sqrt(squareY); 


%WE DO TRANSFORMATIONS FOR HAVING ALIGNED MATRICES:
% Isomorphic Scaling (transformation of a shape smaller or larger while 
% maintaining the ratio of the shapes proportions)
XC = XC / normX;    
YC = YC / normY;

% Rotation
R = XC' * YC;   
[L, D, M] = svd(R);
T = M * L';
    
%Reflection    
M(:,end) = -M(:,end);   
D(end,end) = -D(end,end);
T = M * L';

% Minimized unstandardized distance 
trace = sum(diag(D));   

%Scaling of Y
b = trace * normX / normY;    

%Outputs of the function    
d = 1 - trace.^2;     
Z = normX*trace * YC * T + repmat(columnMeanX, nX, 1);

end 
