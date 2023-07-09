function eta_space = space_estimator(Nodes,Elements,Dirichlet_Edges,u_new,u_old,t, tau, f_h)
% determining tables for edges
[s4e,~,n4s,s4Db,~,e4s] = sides(Elements,Dirichlet_Edges,[]);
% array sizes
nE = size(Elements,1); 
nS = size(n4s,1);
% preallocation as zeros
%eta_S = zeros(nS,1);
eta_T_sq = zeros(nE,1);
eta_3_2_sq = zeros(nE,1);
eta_4_2_sq = zeros(nE,1);
%%Calculate necessary quantities
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
    
b = 1/tau*(u_new-u_old)-f_h; %Used for I.1

%Square of nodal values, pointwise maximum for III.2
u_sq_pntwise_max = max(u_old.^2,u_new.^2);

% loop over all elements
for j = 1:nE
    % element transofrmation 
    X_T = Nodes(Elements(j,:),:); 
    x1 = X_T(1,:);  x2 = X_T(2,:);  x3 = X_T(3,:);
    % Transformation matrix from E_0 to E
    L = [x2-x1; x3-x1; cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1))];
    C = inv(L);
    % Determinant for integral transformation
    dtr = abs(det(L));
    
    % Local mass and stiffness matrix on element E
    M_loc = dtr * M_loc_0;
    
    vol_T = norm(cross(x2-x1,x3-x1))/2; %Area of triangle
    %Formula for Cirmuradius Ris area = product of sides /4R thus the diameter = 2*R
    h_T = norm(x1-x2)*norm(x2-x3)*norm(x3-x1)/(2*vol_T); %Local triangle size

    % I.1 (||1/tau (u_j -u_j-1) -f_j||
    beta = b(Elements(j,:));
    eta_T_sq(j) = h_T^2*beta'*M_loc*beta;
    
    % Term III.2 (squared)
    % The term h_T^2||nalba_Gamma_h u_h ^j|| is pointwise less then maximum of ||nalba_Gamma_h u_h ^j|| and ||nalba_Gamma_h I_j u_h ^j-1||
    A_loc = dtr * ((C(1,1)^2 + C(2,1)^2 + C(3,1)^2)*A_loc11 + (C(1,2)^2 + C(2,2)^2 + C(3,2)^2)*A_loc22 + (C(1,1)*C(1,2) + C(2,1)*C(2,2) + C(3,1)*C(3,2))*(A_loc12+A_loc21));
    beta_pntwise_max = u_sq_pntwise_max(Elements(j,:));
    eta_3_2_sq(j) = h_T^4*beta_pntwise_max'*A_loc*beta_pntwise_max; 
    
    % Term IV.2 (squared)
    beta_f = f_h(Elements(j,:));
    eta_4_2_sq(j) = h_T^4*beta_f'*M_loc*beta_f;

end
eta_S_sq = zeros(nS,1);
eta_S_sq_elementwise = zeros(nE,1);
%% Determine the jump term
% We iterate over all edges instead of all triangles
for j = 1:nS
    % Elements which contain the j-th edge
    T = e4s(j,:);
    if sum(T==0) >0
        continue
    end
    %Define T_1 and T_2
    T_1 = Elements(T(1),:); T_2 = Elements(T(2),:);
    %Permutate the vector such that the third entry is always the node opposite to the shared edge
    while T_1(3) == n4s(j,1) || T_1(3)==n4s(j,2)
        T_1 = circshift(T_1, [0 1]);
        %T_1 = [T_1(end), T_1(1:end-1)];
    end
    while T_2(3) == n4s(j,1) || T_2(3)==n4s(j,2)
        T_2 = circshift(T_2, [0 1]);
        %T_2 = [T_2(end), T_2(1:end-1)];
    end
    % Coordinates of the three points in each trinalge
    T_1_nodes = Nodes(T_1,:); T_2_nodes = Nodes(T_2,:);
    % Compute normals of T_1 and T_2
    AB1 = T_1_nodes(2,:)-T_1_nodes(1,:); BC1 = T_1_nodes(3,:)-T_1_nodes(2,:); CA1 = T_1_nodes(1,:)-T_1_nodes(3,:); %T_1 = conv{ABC}
    BA2 = T_2_nodes(2,:)-T_2_nodes(1,:); AD2 = T_2_nodes(3,:)-T_2_nodes(2,:); DB2 = T_2_nodes(1,:)-T_2_nodes(3,:); %T_2 = conv{BAD} 
    n_T_1 = cross(CA1,-BC1); n_T_1 = n_T_1/norm(n_T_1);
    n_T_2 = cross(DB2,-AD2); n_T_2 = n_T_2/norm(n_T_2);
    % Compute outward Edge normals (n_1_SC refers to the 1-st triangle with edge S
    % opposite to the corner C of the triangle T_1 = conv{A, B, C}
    n_1_SC = -cross(n_T_1,AB1); n_1_SA = -cross(n_T_1,BC1); n_1_SB = -cross(n_T_1,CA1);
    n_2_SD = -cross(n_T_2,BA2); n_2_SB = -cross(n_T_2,AD2); n_2_SA = -cross(n_T_2,DB2);
    %turn to unit normal
    n_1_SC = n_1_SC/norm(n_1_SC); n_1_SB = n_1_SB/norm(n_1_SB); n_1_SA = n_1_SA/norm(n_1_SA); n_2_SD = n_2_SD/norm(n_2_SD); n_2_SB = n_2_SB/norm(n_2_SB); n_2_SA = n_2_SA/norm(n_2_SA);
    % Projection matrices
    P_h_1 = eye(3)-n_T_1' * n_T_1; P_h_2 = eye(3)-n_T_2' * n_T_2;
    % nodal values of the nodes
    u_T_1 = u_new(T_1); u_T_2 = u_new(T_2); %u_T_1 = (u(A) ; u(B); u(C)) u_T_2 = (u(B), u(A), u(D))
    % edge and triangle sizes
    vol_T_1 = norm(cross(AB1,-CA1))/2; vol_T_2 = norm(cross(BA2,-DB2))/2;
    vol_SC_1 = norm(AB1); vol_SA_1 = norm(BC1); vol_SB_1 = norm(CA1);
    vol_SD_2 = norm(BA2); vol_SA_2 = norm(DB2); vol_SB_2 = norm(AD2);
    %sum term over all nodes of the triangle
    sum_1 = -1/(2*vol_T_1) *(u_T_1(1)*vol_SA_1*n_1_SA+ u_T_1(2)*vol_SB_1*n_1_SB+ u_T_1(3)*vol_SC_1*n_1_SC);
    sum_2 = -1/(2*vol_T_2) *(u_T_2(1)*vol_SB_2*n_2_SB +u_T_2(2)*vol_SA_2*n_2_SA+ u_T_2(3)*vol_SD_2*n_2_SD);

    jump = (P_h_1*sum_1')' * n_1_SC' + (P_h_2*sum_2')' * n_2_SD';
    eta_S_sq(j) = jump^2*vol_SC_1^2;
    % Aufteilen der Kantenanteile auf T_1 und T_2 damit wir ein elementwise error haben
    eta_S_sq_elementwise(T) = eta_S_sq_elementwise(T)+eta_S_sq(j);
        
    %Plot triangles with normals
%     CoordsABCD = [T_1_nodes; T_2_nodes(3,:)];
%     trisurf([1 2 3; 2 1 4],CoordsABCD(:,1), CoordsABCD(:,2),CoordsABCD(:,3)) %'FaceColor',['b'])
%     m_AB = T_1_nodes(1,:)+0.5*AB1;
%     m_BC = T_1_nodes(2,:)+0.5*BC1;
%     m_CA = T_1_nodes(3,:)+0.5*CA1;
%     m_AD = T_2_nodes(2,:)+0.5*AD2;
%     m_BA = T_2_nodes(1,:)+0.5*BA2;
%     m_DB = T_2_nodes(3,:)+0.5*DB2;
%     midpoints = [m_AB; m_BC; m_CA; m_BA; m_AD; m_DB];
%     normals = 0.1.*[n_1_SC;n_1_SA;n_1_SB;n_2_SD;n_2_SB;n_2_SA];
%     mid_tri = [1/3*(T_1_nodes(1,:)+T_1_nodes(2,:)+T_1_nodes(3,:)) ;1/3*(T_2_nodes(1,:)+T_2_nodes(2,:)+T_2_nodes(3,:))];
%     mid_norm = -[n_T_1 ; n_T_2];
%     hold on;
%     quiver3(midpoints(1:3,1),midpoints(1:3,2), midpoints(1:3,3), normals(1:3,1), normals(1:3,2), normals(1:3,3),'b')
%     quiver3(midpoints(4:6,1),midpoints(4:6,2), midpoints(4:6,3), normals(4:6,1), normals(4:6,2), normals(4:6,3),'r')
%     quiver3(mid_tri(1:2,1),mid_tri(1:2,2), mid_tri(1:2,3), mid_norm(1:2,1), mid_norm(1:2,2), mid_norm(1:2,3),'k')
%     hold off;
end
eta_space = eta_T_sq+eta_S_sq_elementwise+eta_3_2_sq+eta_4_2_sq; 