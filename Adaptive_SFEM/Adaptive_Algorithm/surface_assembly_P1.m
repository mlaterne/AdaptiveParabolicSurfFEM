function [M,A] = surface_assembly_P1(Nodes,Elements)
% Assembly of mass and stiffness matrix for two dimensional linear finite
% elements on a surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% local degrees of freedom, and full one vector
loc_dof = 3; % due to linear elements
e_loc_dof = ones(1,loc_dof);
% degree of freedom
dof = length(Nodes);


% Insert precomputed local matrices
M_loc_0 = [1/12 1/24 1/24;
        1/24 1/12 1/24;
        1/24 1/24 1/12];
A_loc11 = [1/2 -1/2 0;
           -1/2 1/2 0;
           0 0 0];
A_loc22 = [1/2 0 -1/2;
           0 0 0;
           -1/2 0 1/2];
A_loc12 = [1/2 0 -1/2;
           -1/2 0 1/2;
           0 0 0];
A_loc21 = A_loc12';


%% assembly

M_Cell = cell(size(Elements,1),1);
A_Cell = cell(size(Elements,1),1);

% Assembly loop over elements (elementwise computation)
for i = 1:size(Elements,1)
    
    % element array rows
    Element_array_rows = reshape(e_loc_dof .* Elements(i,:).',[9 1]);
    Element_array_columns = reshape(Elements(i,:) .* e_loc_dof.',[9 1]);
    
    % Load nodes of current elements (3D)
    x1 = [Nodes(Elements(i,1),1) Nodes(Elements(i,1),2) Nodes(Elements(i,1),3)];
    x2 = [Nodes(Elements(i,2),1) Nodes(Elements(i,2),2) Nodes(Elements(i,2),3)];
    x3 = [Nodes(Elements(i,3),1) Nodes(Elements(i,3),2) Nodes(Elements(i,3),3)];
    
    % Transformation matrix from E_0 to E
    L = [x2-x1; x3-x1; cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1))];
    % Matrix of the inverse map
    C = inv(L);
 
    % Determinant for integral transformation
    dtr = abs(det(L));
    
    % Local mass and stiffness matrix on element E
    M_loc = dtr * M_loc_0;
    A_loc = dtr * ((C(1,1)^2 + C(2,1)^2 + C(3,1)^2)*A_loc11 + (C(1,2)^2 + C(2,2)^2 + C(3,2)^2)*A_loc22 + (C(1,1)*C(1,2) + C(2,1)*C(2,2) + C(3,1)*C(3,2))*(A_loc12+A_loc21));
    
    % the vectorized local matrix
    MLoc_vec = reshape(M_loc,[9 1]);
    
    % the vectorized local matrix
    ALoc_vec = reshape(A_loc,[9 1]);
    
    M_Cell{i} = [Element_array_rows Element_array_columns  MLoc_vec];
    A_Cell{i} = [Element_array_rows Element_array_columns ALoc_vec];
    
end

%Cell to sparse matrix
M_Cell_to_sparse = cell2mat( M_Cell );
M = sparse(M_Cell_to_sparse(:,1),M_Cell_to_sparse(:,2),M_Cell_to_sparse(:,3),dof,dof);

A_Cell_to_sparse = cell2mat( A_Cell );
A = sparse(A_Cell_to_sparse(:,1),A_Cell_to_sparse(:,2),A_Cell_to_sparse(:,3),dof,dof);


