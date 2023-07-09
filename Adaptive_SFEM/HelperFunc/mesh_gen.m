function [Nodes, Elements, Dirichlet_Edges,initial_error] = mesh_gen(TOL, fd, boundingbox)
h_0 = TOL; %Needs to be lowered until a Tolerance condition holds

% Surface mode with distmesh to avoid error messages
warning('off','MATLAB:cameramenu:Removal')

Dirichlet_Edges = []; %for closed surfaces

eta_0_global = inf;
while eta_0_global > TOL
    [Nodes,Elements]=distmeshsurface(fd,@huniform,h_0,boundingbox);
    h_0 = h_0 /2;
    eta_0_global = h_0;
end

initial_error = eta_0_global;

disp('Initialsation done')
