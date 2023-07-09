%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Problem
% This script contains the initialization of the numerical Experiments
%
% T_0 = Start_time
% T_end = Final time
% theta_refine = in [0, 1]; %Used for Refinment -> Large Theta = Few refinements
% theta_coarse = in [0, 1]; %Used for Coarsening -> Large Theta = Many coarsenings
% TOL = Squared Tolerance bounding the error
% h_0 = Initial Mesh Size (used by distmesh)
% Nodes_0 = Initial Node Set
% Elements_0 = Initial set of Elements
% Dirichlet_Edges_0 = Initial Dirichlet Edges (in all upcoming experiment = \emptyset)
% Type of Refinement: [1=RGB 2=NVB] 
% Type of Coarsen: [No Coarsening=0 RGB=1 NVB=2 Naiv=3]
% Marking Criterion: {1 = Bulk, 2 = Dörfler]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set All Parameters:
T_0 = 0; T_end = 1; % Time Interval
coarsen_type = 1; refinement_type = 1; marking_criterion = 1; % Marking, Refinement, Coarsening Strategy
theta_refine = 0.5; theta_coarse = 0.25; %Refinement Parameters 
h_0 = 0.8; %Initial Elementsize used by distmesh
TOL = 1; %Tol for adaptive algorithm

%% Initialize Surface
d = @(x) -(1- (x(:,1).^2)-x(:,2).^2-x(:,3).^2); %Oriented Distance function for sphere
bb = 1.2*[-1,-1,-1;1,1,1]; %Bounding Box of the Surface

%% PDE
f = @(x,t) 5* exp(-t) * x(:,1).* x(:,2); % the functions u = x_j x_i (i \neq j = 1,2,3) are eigenfunction of the Laplace-Beltrami operator on the unit sphere, eigenvalue (-6)
u = @(x,t) exp(-t) * x(:,1).* x(:,2); %Solution if known (only needed for comparison)
u_0 = @(x) x(:,1).* x(:,2); %Initial condition

% Mesh 0
[Nodes_0, Elements_0] = mesh_gen(h_0,d,bb); %Requires a surface parametrasized by distance function d and bounding box 
Dirichlet_Edges_0 = []; %Surface is closed
Nodes_0 = lift(Nodes_0,d); %Higher Accuracy lift

%% Run Code
OUTPUT = Main(TOL, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d ,f, u_0,coarsen_type, refinement_type,marking_criterion);
%%%%%%%%%%%%%%
% The Output is a NxM cell which contains for the discrete steps j = 1,...,N multiple values
% The OUTPUT coloumns are:
% 1. Nodes
% 2. Elements
% 3. Dirichlet_Edges
% 4. u_h ^j Nodal values of approximation for j-th timestep
% 5. u_h ^{j-1 -> j} the solution for the previous timestamp interpolated onto the new mesh
% 6. Time t_j
% 7. Global Error indicated by spatial ndicator
% 8. Global Error indicated by temporal indicator
% 9. Right-hand-side f interpolated on \Gamma_j
% 10. Mass matrix
% 11. Stiffness Matrix
%%%%%%%%%%%%%%

%Plot the Surface at the final timestep, colored according to u_h ^j
Plot(cell2mat(OUTPUT(end,4)),cell2mat(OUTPUT(end,2)),cell2mat(OUTPUT(end,1)))



function hh = Plot(func, Elements, Nodes)
%Plot Surface and Color faces according to func 
C_Tri = func; 
hh = trisurf(Elements,Nodes(:,1),Nodes(:,2),Nodes(:,3));
 % Additional bit to control color of each patch individually
set(gca,'CLim',[min(C_Tri), max(C_Tri)]);
set(hh,'FaceColor','flat',...
'FaceVertexCData',C_Tri,...
'CDataMapping','scaled');
axis equal;
xlabel('x') 
ylabel('y') 
zlabel('z')
end


