function u_new = solve_SFEM(u_old, Nodes, Dirichlet_Edges, M, A, time, tau, f_h)
innerNodes = setdiff(1:size(Nodes,1),unique(Dirichlet_Edges)); %Alle Knoten die nicht am Rand sind

%Solve (M+\tau A) u_j = M (u_(j-1)+ \tau *f_func_j)
lhs = (M + tau * A);

%t_new = time + tau;
    
% r.h.s.
b = f_h; %f_func_j
rhs = M * (tau * b + u_old);
    
% solving the linear system
u_new = zeros(size(Nodes,1),1);
u_new(innerNodes) = lhs(innerNodes,innerNodes) \ rhs(innerNodes);


end