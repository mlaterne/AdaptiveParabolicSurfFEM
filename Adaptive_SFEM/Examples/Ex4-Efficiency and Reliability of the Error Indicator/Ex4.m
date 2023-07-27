%% Example 4: Ohne Adaptivität für ein Satz von Gittern und Zeitschrittweiten L2 Estimator und L2 Error plotten

%% Initialization from Ex1
T_0 = 0; T_end = 1;
%% Initialize Surface
d = @(x) -(1- (x(:,1).^2)-x(:,2).^2-x(:,3).^2); %Oriented Distance function for sphere
Dirichlet_Edges = []; neumann = [];
s = sqrt(2)/2;
Nodes_0 = [0 0 -1; -s -s 0; s -s 0; s s 0; -s s 0; 0 0 1]; %Ocathedron
Elements_0 = [1 2 3; 1 3 4; 2 1 5; 5 1 4; 3 2 6; 4 3 6; 2 5 6; 5 4 6];
%%PDE
f = @(x,t) 5* exp(-t) * x(:,1).* x(:,2); % the functions u = x_j x_i (i \neq j = 1,2,3) are eigenfunction of the Laplace-Beltrami operator on the unit sphere, eigenvalue (-6)
u = @(x,t) exp(-t) * x(:,1).* x(:,2); %Solution
u_0 = @(x,t) x(:,1).* x(:,2); %Initial condition

%% Set of Time steps and meshes
tau_all = [1 0.1 0.01];
mesh = cell(2,7); mesh(:,1) = {Nodes_0, Elements_0}; DOFs(1) = size(Nodes_0,1);
for i = 1:6
    marked = logical(ones(size(cell2mat(mesh(2,i)),1),1));
    [Nodes,Elements,~, Dirichlet_Edges] = TrefineRGB(cell2mat(mesh(1,i)),cell2mat(mesh(2,i)),Dirichlet_Edges,neumann,marked);
    Nodes = lift(Nodes,d);
    mesh(:,i+1) = {Nodes, Elements};
    DOFs(i+1) = size(Nodes,1);
end

%Compute the nonadaptive solution and extract error indicators and l2 and
%h1 error quantaties for all pairs of (tau,mesh)
for j = 1:3
    tau = tau_all(j);
    for i = 1:7
        OUTPUT = NoAdaptivitySFEMPara(cell2mat(mesh(1,i)), cell2mat(mesh(2,i)), Dirichlet_Edges, T_0, T_end ,f, u_0,tau);
        eta_space = cell2mat(OUTPUT(:,7)); eta_temp = cell2mat(OUTPUT(:,8));
        L2_Estimator(i,j) = max(eta_space+eta_temp);
        [L2_Error(i,j), H1_Error(i,j)] = L2_Ex1(OUTPUT,u);
        display([string(datetime)+ "     Solution found for the index pair (i,j) = (" + num2str(i)+ ","+ num2str(j)+ ")"])
    end
end

%% Plot Error curves
f5 = figure;
symbols='soxd+^*v><';
subplot(1,2,1)
for j = 1:3
    loglog(DOFs,L2_Error(1:end,j),'Marker',symbols(j),'LineStyle','-','color','k');
    hold on
    loglog(DOFs,1.*L2_Estimator(1:end,j),'Marker',symbols(j),'LineStyle','-','color','b');
end
xlabel("DOF")
loglog(DOFs,4./DOFs,'r--')
loglog(DOFs,1200./DOFs,'r--')
legend('L2 Error $(\tau = 1)$', 'Estimator $(\tau = 1)$','L2 Error $(\tau = 0.1)$', 'Estimator $(\tau = 0.1)$','L2 Error $(\tau = 0.01)$', 'Estimator $(\tau = 0.01)$', '$\mathcal{O} (DOF^{-1})$','Interpreter', 'latex','FontSize', 8)
title('Estimator and L2 Error-term for different $\tau$','Interpreter', 'latex')
%hold off
subplot(1,2,2)
for j = 1:3
    %loglog(DOFs,L2_Error(1:end,j),'Marker',symbols(j),'LineStyle','-','color',[0.3 0.3 0.3]);
    loglog(DOFs,H1_Error(1:end,j),'Marker',symbols(j),'LineStyle','-','color','k');
    hold on
    loglog(DOFs,1.*L2_Estimator(1:end,j),'Marker',symbols(j),'LineStyle','-','color','b');
end
xlabel("DOF")
loglog(DOFs,14./DOFs,'r--')
legend('H1 Error $(\tau = 1)$', 'Estimator $(\tau = 1)$','H1 Error $(\tau = 0.1)$', 'Estimator $(\tau = 0.1)$','H1 Error $(\tau = 0.01)$', 'Estimator $(\tau = 0.01)$', '$\mathcal{O} (DOF^{-1})$','Interpreter', 'latex','FontSize', 8)
title('Estimator and H1 Error-term for different $\tau$','Interpreter', 'latex')
hold off
set(gcf, 'Position', [100, 100, 1000, 500]);

%% Effect of Lift
surfmesh_quality_ex4(Nodes, Elements)
        
 function [L2, H1] = L2_Ex1(OUTPUT,u) 
%Plot L2 Error
for j = 1:size(OUTPUT,1)
    u_h = cell2mat(OUTPUT(j,4));
    Nodes = cell2mat(OUTPUT(j,1));
    t = cell2mat(OUTPUT(j,6));
    %% computing and storing errors
    if j == 1
        errM(j) = 0;
        errH(j) = 0;
    else
        err = u_h - u(Nodes,t);
        M = OUTPUT{j,10}; A = OUTPUT{j,11};
        errM(j) = sqrt(err' * M * err);
        errH(j) = err' * A * err;
    end
end
L2 = max(errM);
t_all = cell2mat(OUTPUT(:,6));
H1 = sqrt(sum(0.5.*(errH(2:end)+errH(1:(end-1)))'.*(t_all(2:end)-t_all(1:(end-1)))));
 end

function surfmesh_quality_ex4(Nodes, Elements)

%% mesh quality
for ind=1:length(Elements)
    % the three nodes
    A=Nodes(Elements(ind,1),:);
    B=Nodes(Elements(ind,2),:);
    C=Nodes(Elements(ind,3),:);
    
    % ratio of edge lengths
    ratio_Nodes(ind)=max([norm(A-B) norm(C-B) norm(A-C) ])/min([norm(A-B) norm(C-B) norm(A-C) ]);
    
    % angles
    angles=[ acosd(dot((B-A)/norm(B-A),(C-A)/norm(C-A)))  acosd(dot((A-B)/norm(A-B),(C-B)/norm(C-B)))  acosd(dot((A-C)/norm(A-C),(B-C)/norm(B-C))) ];
    angle_min_Nodes(ind)=min(angles);
    angle_max_Nodes(ind)=max(angles);
    
    % normalized angle ratio
    skewness_Nodes(ind)=max( (angle_max_Nodes(ind)-60)/120 , (60-angle_min_Nodes(ind))/60 );
    
    % inscribed and outscribed circle ratio
    s=0.5*( norm(A-B)+norm(A-C)+norm(B-C) );
    ic=sqrt( s*(s-norm(A-B))*(s-norm(A-C))*(s-norm(B-C)) )/s;
    cc=norm(A-B)/( 2*sin(acos(dot((A-C)/norm(A-C),(B-C)/norm(B-C)))) );
    circ_ratio_Nodes(ind)=cc/ic;
    
end

%% plotting mesh quality measures
fctr=1;

f2 = figure;
subplot(2,2,1)
histogram(angle_min_Nodes)
ylabel('# of Elements')
title('Minimum angles')


subplot(2,2,3)
histogram(angle_max_Nodes)
ylabel('# of Elements')
title('Maximum angles')

subplot(2,2,2)
histogram(skewness_Nodes)

ylabel('# of Elements')
title('Normalized Angle Ratio (the lower the better)')

subplot(2,2,4)
histogram(circ_ratio_Nodes)
ylabel('# of Elements')
title('Inscribed and Circumscribed Circle Radius Ratio')
end