function OUTPUT = NoAdaptivitySFEMPara(Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end ,f, u_0,tau)

%% Solving Parabolic Surface PDEs with Adaptiv FEM (based on residual L2 estimator)
% This is the main script which calls all kinds of functions (calculating
% the approximated and exact solution, plotting of residual and full L2
% error and many more.

% The script is divided into multiple steps, based on Bartels adaptivity
% algorithm. Some functions were provided by Kovács, Bartels and Funken

%% Initializing data

%Initialise problem, and store inital conditons into solution space
t = T_0;
j= 1;
f_h = f(Nodes_0,T_0);
u_old = u_0(Nodes_0);

OUTPUT(1,:) = {Nodes_0,Elements_0,Dirichlet_Edges_0,u_old,u_old,T_0,0, 0,f_h,zeros(size(Nodes_0,1)),zeros(size(Nodes_0,1))};

%Nodes_0 store inital grid
Nodes = Nodes_0; Elements = Elements_0;Dirichlet_Edges = Dirichlet_Edges_0;

while t < T_end
    t = t + tau; %New timestep
    
    [M,A] = surface_assembly_P1(Nodes,Elements); % Mass and stiffness matrices
    % The implicit scheme gives (M+\tau A) u_j = M (u_(j-1)+ \tau *f_func_j)
    f_h = f(Nodes,t);
    u_new = solve_SFEM(u_old, Nodes, Dirichlet_Edges, M, A, t, tau, f_h);
        
    %% Determine spatial error component
    eta_space_sq = space_estimator(Nodes,Elements,Dirichlet_Edges,u_new,u_old,t, tau, f_h); %Berechnet den Lokalen Fehlerindikator (jeder Eintrag steht für den Fehler auf dem Element)
    eta_space_global = sum(eta_space_sq); %Gesamtfehler

%% Compute Temporal Error
    eta_temp_sq = temp_estimator(Nodes,Elements,Dirichlet_Edges,u_new, u_old ,t,tau,f);
    eta_temp_global = sum(eta_temp_sq);

    if j == 1 %Overwrite the initial mesh with the one refined in the first step
        OUTPUT(1,:) = {Nodes,Elements,Dirichlet_Edges,u_old,u_old,T_0,0, 0,f_h,[],[]};
    end
    if t+tau > T_end
        tau = T_end-t;
    end
    j = j+1;
    %Speicher Lösung und Quantities zum Plotten in einem Cell array
    OUTPUT(j,:) = {Nodes,Elements,Dirichlet_Edges,u_new, u_old ,t,eta_space_global, eta_temp_global, f_h, M, A}; %Also store u_old which is already nodally interpolated to be able to write the solutions for continous t values 
        
    %Update das die neue Lösung jetzt die alte ist
    u_old = u_new;
        
    %disp(['Part of solution found at t = ', num2str(t)])     
end
