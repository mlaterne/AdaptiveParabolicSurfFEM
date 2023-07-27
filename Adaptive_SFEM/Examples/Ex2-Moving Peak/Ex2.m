%% Experiment 2: Moving Peak
T_0 = 0; T_end = 1; theta_refine = 0.5; theta_coarse = 0.25; set_TOL = 0.18;
coarsen_type = 2; refinement_type = 2; marking_criterion = 1;
% Initialize Surface
d = @(x) -(1- (x(:,1).^2)-x(:,2).^2-x(:,3).^2); %Oriented Distance function for sphere
grad_d = @(x) 2.*[x(:,1), x(:,2) ,x(:,3)];
bb = 1.2*[-1,-1,-1;1,1,1]; %Bounding Box of the Surface
[Nodes_0, Elements_0, Dirichlet_Edges_0, initial_error] = mesh_gen(set_TOL,d, bb); %Requires a surface parametrasized by distance function d and bounding box 
Nodes_0 = lift(Nodes_0,d);
% PDE
a = 25; b = 600; T = 2;

u = @(x,t) (1-exp(-b*(t-0.5).^2)).*exp(-a*((x(:,1)-cos(t./T *pi)).^2 + (x(:,2)-sin(t./T .*pi)).^2 + x(:,3).^2));
u_0 = @(x) u(x,0);
f = @(x,t) -8.*((1./2.*T.*a.^2.*(x(:,1)-x(:,2)).*(x(:,1)+x(:,2)).*cos(t.*pi./T).^2+a.*(sin(t.*pi./T).*T.*a.*x(:,1).*x(:,2)+1./2.*T.*x(:,1)+1./4.*pi.*x(:,2)).*cos(t.*pi./T)+1./2.*(x(:,2).*T-1./2.*pi.*x(:,1)).*a.*sin(t.*pi./T)+1./2.*T.*((x(:,2).^2-1).*a.^2-1./2.*(t-1./2).*b)).*exp(-1./4.*b+b.*t)-(1./2.*T.*a.*(x(:,1)-x(:,2)).*(x(:,1)+x(:,2)).*cos(t.*pi./T).^2+(sin(t.*pi./T).*T.*a.*x(:,1).*x(:,2)+1./2.*T.*x(:,1)+1./4.*pi.*x(:,2)).*cos(t.*pi./T)+(1./2.*x(:,2).*T-1./4.*pi.*x(:,1)).*sin(t.*pi./T)+1./2.*T.*a.*(x(:,2)-1).*(x(:,2)+1)).*a.*exp(b.*t.^2)).*exp(2.*a.*x(:,2).*sin(t.*pi./T)+2.*a.*x(:,1).*cos(t.*pi./T)-b.*t.^2-2.*a)./T;

%% Compute approximation and plot the L2 error
OUTPUT = Main(1, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d , f , u_0,coarsen_type, refinement_type,marking_criterion);
%surfmesh_quality(OUTPUT{1,1}, OUTPUT{1,2})
[a,L2] = L2_Ex1(OUTPUT,u);

%% Plot the Peak for different times
steps = [1, 127, 142, 285];
fig = figure;
i = 1;
for j = steps
    subaxis(2,2,i,'Spacing', 0, 'Padding', 0.02, 'Margin', 0.03)
    i = i+1;
    t = OUTPUT{j,6}; %Midpoint at j = 35;
    Nodes = cell2mat(OUTPUT(j,1)); Elements = cell2mat(OUTPUT(j,2)); u_new = cell2mat(OUTPUT(j,4));
    C_Tri = u_new; 
    hh =trisurf(Elements,Nodes(:,1),Nodes(:,2),Nodes(:,3));
    set(hh, 'FaceColor', 'white');
    axis equal;
    xlabel('x') 
    ylabel('y') 
    zlabel('z')
    title(['t = ' , num2str(t)])
    view(135, 15)
end

% %% Timeit Test 
% a = 5; b = 100; T = 1;
% u = @(x,t) (1-exp(-b*(t-0.5).^2)).*exp(-a*((x(:,1)-cos(t./T *pi)).^2 + (x(:,2)-sin(t./T .*pi)).^2 + x(:,3).^2));
% u_0 = @(x) u(x,0);
% f = @(x,t) -8.*((1./2.*T.*a.^2.*(x(:,1)-x(:,2)).*(x(:,1)+x(:,2)).*cos(t.*pi./T).^2+a.*(sin(t.*pi./T).*T.*a.*x(:,1).*x(:,2)+1./2.*T.*x(:,1)+1./4.*pi.*x(:,2)).*cos(t.*pi./T)+1./2.*(x(:,2).*T-1./2.*pi.*x(:,1)).*a.*sin(t.*pi./T)+1./2.*T.*((x(:,2).^2-1).*a.^2-1./2.*(t-1./2).*b)).*exp(-1./4.*b+b.*t)-(1./2.*T.*a.*(x(:,1)-x(:,2)).*(x(:,1)+x(:,2)).*cos(t.*pi./T).^2+(sin(t.*pi./T).*T.*a.*x(:,1).*x(:,2)+1./2.*T.*x(:,1)+1./4.*pi.*x(:,2)).*cos(t.*pi./T)+(1./2.*x(:,2).*T-1./4.*pi.*x(:,1)).*sin(t.*pi./T)+1./2.*T.*a.*(x(:,2)-1).*(x(:,2)+1)).*a.*exp(b.*t.^2)).*exp(2.*a.*x(:,2).*sin(t.*pi./T)+2.*a.*x(:,1).*cos(t.*pi./T)-b.*t.^2-2.*a)./T;
% stratgies = [1 1; 2 2; 1 3; 2 3; 1 0; 2 0]; %Refining, Coarsening
% for i = 1:1
%     myFunction = @() Main(1, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d, f , u_0,  stratgies(i,2), stratgies(i,1),marking_criterion);
%     execution_time(i) = timeit(@() cellfun(@(x) [], myFunction(), 'UniformOutput', false));
% end


function [a,L2, b, H1] = L2_Ex1(OUTPUT,u) 
%Plot L2 Error
for j = 1:size(OUTPUT,1)
    u_h = cell2mat(OUTPUT(j,4));
    Nodes = cell2mat(OUTPUT(j,1));
    t = cell2mat(OUTPUT(j,6));
    %% computing and storing errors
    err = u_h - u(Nodes,t);
    if j == 1
        errM(j) = 0;
        errH(j) = 0;
    else
        M = OUTPUT{j,10}; A = OUTPUT{j,11};
        errM(j) = sqrt(err.' * M * err);
        errH(j) = err.' * A * err;
        L2 = max(errM);

    end
end
t_all = cell2mat(OUTPUT(:,6));
H1 = sum(0.5.*(errH(2:end)+errH(1:(end-1)))'.*(t_all(2:end)-t_all(1:(end-1))));
a = plot(cell2mat(OUTPUT(:,6)), errM);
b = 1;
end

function hh = Plot(eta_space, Elements, Nodes)
%Plot Elementwise Error 
C_Tri = eta_space; 
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