function eta_temp = temp_estimator(Nodes,Elements,Dirichlet_Edges,u_new, u_old ,t,tau, f)
%% Compute the temporal error indicator, only usable for closed surfaces

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
 
% array sizes
[n,m] = size(Nodes);
nE = size(Elements,1); 

%nodal difference of the current solution and the old interpolated
%solution:
b = u_new-u_old;

% preallocation as zeros eta_x is the x-th term we have to compute
eta_2 = zeros(nE,1); %h^4 /tau^2 ||u_h ^j - I_j u_h ^j-1||^2
eta_3 = zeros(nE,1); %||nabla (u_h ^j - I_j u_h ^j-1)||^2
eta_4 = zeros(nE,1); %Determine tau^2 sup ||d_t f||^2
% loop over all elements
for j = 1:nE
    % element transofrmation 
    X_T = Nodes(Elements(j,:),:); 
    x1 = X_T(1,:);  x2 = X_T(2,:);  x3 = X_T(3,:);
    % Transformation matrix from E_0 to E
    L = [x2-x1; x3-x1; cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1))];
    % Matrix of the inverse map
    C = inv(L);
 
    % Determinant for integral transformation
    dtr = abs(det(L));
    
    % Local mass and stiffness matrix on element E
    M_loc = dtr * M_loc_0;

    vol_T = norm(cross(x2-x1,x3-x1))/2; %Area of triangle
    %Formula for Cirmuradius Ris area = product of sides /4R thus the diameter = 2*R
    h_T = norm(x1-x2)*norm(x2-x3)*norm(x3-x1)/(2*vol_T); %Local triangle size
    
    %% (II)-Term
    % On each triangle we have to calculate ||u_h ^j - I_j u_h^(j-1)||_{L^2(T)}
    % Using the nodal value b_i = u_i^j - I_j u_h^(j-1):
    beta = b(Elements(j,:));
    eta_2(j) = h_T^4 *beta'*M_loc*beta;
    
    %% (III)-Term
    % Determine ||nalba_Gamma_h (u_h ^j - I_j u_h^(j-1))||_{L^2(T)}
    A_loc = dtr * ((C(1,1)^2 + C(2,1)^2 + C(3,1)^2)*A_loc11 + (C(1,2)^2 + C(2,2)^2 + C(3,2)^2)*A_loc22 + (C(1,1)*C(1,2) + C(2,1)*C(2,2) + C(3,1)*C(3,2))*(A_loc12+A_loc21));
    eta_3(j) = beta'*A_loc*beta;
    
    %% (IV)-Term
    % Determine sup ||f(t_j)-f(t)||+h_T^2 (t in (t_j-tau,t_j))
    center = (x1+x2+x3)/3; % The centroid of the triangle
    f_t = 1/(eps^(1/2)) *(f(center,t+eps^(1/2))- f(center,t));
    if abs(f_t) == Inf
        warning(["The function f evaluated to be Inf for (x,t) = " + "("+ num2str(center)+ ","+ num2str(t)+ ")"])
        f_t = 0;
    end
    eta_4(j) = vol_T*f_t^2;
end
eta_temp = 1/tau^2 *eta_2+ eta_3+ tau^2 *eta_4;