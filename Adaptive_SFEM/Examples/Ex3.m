%% Experiment 3: Pertubated Torus - Mesh Quality
T_0 = 0; T_end = 250; theta_refine = 0.5; theta_coarse = 0.25; set_TOL = 0.3;
coarsen_type = 1; refinement_type = 1; marking_criterion = 1;

%% Initialize Surface
R = 3; r = 0.5;
d= @(x) ( (x(:,1).^2 + x(:,2).^2).^(0.5) - R ).^2 + x(:,3).^2 - ( r + 0.25 * sin( 6 * atan2(x(:,2),x(:,1)) ) ).^2; 
grad_d = @(x) 0;
bb = 12*[-1,-1,-0.2;1,1,0.2]; %Bounding Box of the Surface
[Nodes_0, Elements_0, Dirichlet_Edges_0, initial_error] = mesh_gen(set_TOL,d, bb); %Requires a surface parametrasized by distance function d and bounding box 
Nodes_0 = lift(Nodes_0,d);
%%PDE
f = @(x,t) 0 + 0* exp(-t) .* x(:,1).* x(:,2); %Zero function
u_0 = @(x) 1 ./ (1 + exp(-100.*x(:,1)));

% Plot Solution after first timestep
OUTPUT = Main(12, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d ,f, u_0,coarsen_type, refinement_type,marking_criterion);
Nodes_1 = cell2mat(OUTPUT(2,1)); Elements_1 = cell2mat(OUTPUT(2,2)); u_1= cell2mat(OUTPUT(2,4));
Plot(u_1, Elements_1, Nodes_1)

%% Plot the torus for multiple times
% fix the color scale and plot multiple time steps
C_Tri_fix = cell2mat(OUTPUT(1,4));
k = 1;
fig = figure();
for i = [1 7 13 size(OUTPUT,1)]
    Nodes = cell2mat(OUTPUT(i,1)); Elements = cell2mat(OUTPUT(i,2)); C_Tri= cell2mat(OUTPUT(i,4)); t = cell2mat(OUTPUT(i,6));     
    subplot(2,2,k); k = k+1;
    hh = trisurf(Elements,Nodes(:,1),Nodes(:,2),Nodes(:,3)); hh.LineWidth = 0.05;
 % Additional bit to control color of each patch individually
    ax = gca;
    ax.Colormap = colormap(parula(256));
    set(ax,'CLim',[min(C_Tri_fix), max(C_Tri_fix)]);
    set(hh,'FaceColor','flat',...
    'FaceVertexCData',C_Tri,...
    'CDataMapping','scaled');
    axis equal;
    xlabel('x') 
    ylabel('y') 
    zlabel('z')
    title(['t = ' , num2str(t)])
    view(90,90)
end
%Break and add the colorbar manually
h = axes(fig,'visible','off'); 
c = colorbar(h,'Position',[0.9300 0.1050 0.0400 0.8300]);colormap(parula(256)); caxis(h,[min(C_Tri_fix), max(C_Tri_fix)]); %cb.Layout.Tile = 'east';
hold off
fig2 = figure();
Nodes = cell2mat(OUTPUT(1,1)); Elements = cell2mat(OUTPUT(1,2)); C_Tri= cell2mat(OUTPUT(1,4)); t = cell2mat(OUTPUT(1,6));     
hh = trisurf(Elements,Nodes(:,1),Nodes(:,2),Nodes(:,3)); hh.LineWidth = 0.05;

%% Mesh Quality
coarsen_type = 1; refinement_type = 1; marking_strategy = 1;
OUTPUT = Main(12, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d ,f, u_0,coarsen_type, refinement_type,marking_criterion);
Nodes_1 = cell2mat(OUTPUT(2,1)); Elements_1 = cell2mat(OUTPUT(2,2)); u_1= cell2mat(OUTPUT(2,4));
surfmesh_quality(Nodes_0, Elements_0)
surfmesh_quality(Nodes_1, Elements_1)



function hh = Plot(func, Elements, Nodes)
%Plot Elementwise Error 
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


function surfmesh_quality(Nodes, Elements)

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
line([20 20], ylim, 'Color', 'k', 'LineWidth', 2,'LineStyle', '--')
ylabel('# of Elements')
title('Minimum angles')


subplot(2,2,3)
histogram(angle_max_Nodes)
line([120 120], ylim, 'Color', 'k', 'LineWidth', 2,'LineStyle', '--')
ylabel('# of Elements')
title('Maximum angles')

subplot(2,2,2)
histogram(skewness_Nodes)
line([.5 .5], ylim, 'Color', 'k', 'LineWidth', 2,'LineStyle', '--')
line([.8 .8], ylim, 'Color', 'k', 'LineWidth', 2,'LineStyle', ':')
ylabel('# of Elements')
title('Normalized Angle Ratio (the lower the better)')

subplot(2,2,4)
histogram(circ_ratio_Nodes)
line([10 10], ylim, 'Color', 'k', 'LineWidth', 2,'LineStyle', '--')
ylabel('# of Elements')
title('Inscribed and Circumscribed Circle Radius Ratio')
end