%% Experiment 1: Exponentially Decaying Solution with multiple figures
%% Example 1:
T_0 = 0; T_end = 1; theta_refine = 0.5; theta_coarse = 0.25; set_TOL = 0.8;
coarsen_type = 1; refinement_type = 1; marking_criterion = 1;
%% Initialize Surface
d = @(x) -(1- (x(:,1).^2)-x(:,2).^2-x(:,3).^2); %Oriented Distance function for sphere
grad_d = @(x) 2.*[x(:,1), x(:,2) ,x(:,3)];
bb = 1.2*[-1,-1,-1;1,1,1]; %Bounding Box of the Surface
[Nodes_0, Elements_0, Dirichlet_Edges_0, initial_error] = mesh_gen(set_TOL,d, bb); %Requires a surface parametrasized by distance function d and bounding box 
Nodes_0 = lift(Nodes_0,d);
%%PDE
f = @(x,t) 5* exp(-t) * x(:,1).* x(:,2); % the functions u = x_j x_i (i \neq j = 1,2,3) are eigenfunction of the Laplace-Beltrami operator on the unit sphere, eigenvalue (-6)
u = @(x,t) exp(-t) * x(:,1).* x(:,2); %Solution
u_0 = @(x,t) x(:,1).* x(:,2); %Initial condition

%% Plot the L2 Norm for a set of large Tolerances
TOL = [25 22 20 17.5 15 12 9 7 3.5 1]; 
for i = 1:size(TOL,2)
    OUTPUT = Main(TOL(i),Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d, f , u_0, coarsen_type,refinement_type,marking_criterion);
    Z(i) = {OUTPUT};
end
% L2 Error over time
for i = 1:size(TOL,2)
    [~,L2(i),~,H1(i)] = L2_Ex1(Z{i},u);
end 

% L2 Error against TOL
f_L2 = figure();
hold on
plot([0 TOL(1)], [0 L2(1)],'k--')
x_lin = linspace(0,0.5,1000); x_sqr = sqrt(x_lin)*0.009; 
%plot(x_lin, x_sqr,'b--')
scatter(TOL, L2,50,'x','r')
legend({'$\mathcal{O} (TOL)$', '$\approx \sup_t \| e(t)\|$'},'Interpreter','latex','Location','northwest','FontSize', 12)
xlabel('TOL')
ylabel('Sup of L2-Error')
%title(['Supremum of L2 Error vs. TOL'])
hold off

% H1 Error against TOL
f_H1 = figure();
hold on
%plot([0 TOL(1)], [0 H1(1)],'k--')
x_lin = linspace(0,TOL(1),1000); 
plot(x_lin, 0.003+0.000143*x_lin.^2,'b--')
scatter(TOL, H1,50,'x','r')
legend({'$\mathcal{O} (TOL ^2)$', '$\approx \int_0 ^T \| \nabla_{\Gamma} e(t) \|_{L^2} ^2$'},'Interpreter','latex','Location','northwest','FontSize', 12)
xlabel('TOL')
ylabel('Integral of H1-Seminorm Error squared')
%title(['Supremum of L2 Error vs. TOL'])
hold off

%% Plot the L2 Norm for a set of small Tolerances (takes some time)
% TOL = [0.5 0.4 0.3 0.2 0.1 0.05 0.035 0.02];
% clear L2, clear H1
% % for i = 1:size(TOL,2)
% %     OUTPUT = Main(TOL(i),Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d, f , u_0, coarsen_type,refinement_type,marking_criterion);
% %     Z(i) = {OUTPUT};
% % end
% 
% for i = 1:size(TOL,2)
%     [~,L2(i),~,H1(i)] = L2_Ex1(Z{i},u);
% end 
% 
% f2_L2 = figure();
% hold on
% x_lin = linspace(0,0.5,1000); x_sqr = sqrt(x_lin)*0.009; 
% plot(x_lin, x_sqr,'k--')
% scatter(TOL, L2,50,'x','r')
% legend({'$\mathcal{O} (\sqrt{TOL})$', '$\approx \sup_t \| e(t)\|$'},'Interpreter','latex','Location','northwest','FontSize', 12)
% xlabel('TOL')
% ylabel('Sup of L2-Error')
% hold off
% 
% f2_H1 = figure();
% hold on
% x_lin = linspace(0,0.5,1000); x_sqr = sqrt(x_lin)*0.002; 
% plot(x_lin, x_lin*0.005,'b--')
% scatter(TOL, H1,50,'x','r')
% legend({'$\mathcal{O} (TOL)$', '$\approx \int_0 ^T \| \nabla_{\Gamma} e(t) \|_{L^2} ^2$'},'Interpreter','latex','Location','northwest','FontSize', 12)
% xlabel('TOL')
% ylabel('Integral of H1-Seminorm Error squared')
% hold off

%% Plot Nodes vs Time
TOL = 0.1; 
T_end = 3;
OUTPUT = Main(TOL,Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d, f , u_0, coarsen_type,refinement_type,marking_criterion);

Nodes = OUTPUT(:,1);
t = cell2mat(OUTPUT(:,6));
for i = 1:size(Nodes,1)
    numb_nodes(i) = size(Nodes{i},1);
end
figure;
stairs(t,numb_nodes,'k','LineWidth',3)
ylabel('# of Nodes')
xlabel('t')
set(gca,'YScale','log');

%Plot L2 Error over time
figure;
[t,errM] = L2_Ex1_Plot(OUTPUT,u);
plot(t,errM,'k')
title({'Approximation to the pointwise $L^2$ error by $\|u_h^j - \tilde{I_h} u^\ell (t_j) \|_{L^2 (\Gamma_h)} = \sqrt{(\mathbf{e}^j) ^T \mathbf{M} \mathbf{e}^j}$'},'Interpreter','latex')
xlabel('t')

%eta_space and eta_temp over time
figure;
eta_space = cell2mat(OUTPUT(2:end,7));
eta_temp = cell2mat(OUTPUT(2:end,8));
plot(t(2:end), eta_space, 'LineWidth', 1.5)
hold on
plot(t(2:end), eta_temp, 'LineWidth', 1.5)
legend('\eta_{space}', '\eta_{temp}','Intepreter','latex', 'Location', 'northeast','FontSize', 12)
xlabel('t','FontSize', 12)
ylabel('Indicated Error','FontSize', 12)
hold off


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
    end
end
L2 = max(errM);
t_all = cell2mat(OUTPUT(:,6));
H1 = sum(0.5.*(errH(2:end)+errH(1:(end-1)))'.*(t_all(2:end)-t_all(1:(end-1))));
a = 1;
b = 1;
end

function [t,errM] = L2_Ex1_Plot(OUTPUT,u) 
%Plot L2 Error
for j = 1:size(OUTPUT,1)
    u_h = cell2mat(OUTPUT(j,4));
    Nodes = cell2mat(OUTPUT(j,1));
    t = cell2mat(OUTPUT(j,6));
    %% computing and storing errors
    err = u_h - u(Nodes,t);
    if j == 1
        errM(j) = 0;
    else
        M = OUTPUT{j,10};
        errM(j) = sqrt(err.' * M * err);
    end
end
t = cell2mat(OUTPUT(:,6));
end