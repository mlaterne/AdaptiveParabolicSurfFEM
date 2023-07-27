function OUTPUT = Main(set_TOL, Nodes_0, Elements_0, Dirichlet_Edges_0, T_0, T_end, theta_refine, theta_coarse, d ,f, u_0,coarsen_type, refinement_type,marking_criterion)

%% Solving Parabolic Surface PDEs with Adaptiv FEM (based on residual L2 estimator)
% This is the main script which calls all kinds of functions (calculating
% the approximated and exact solution, plotting of residual and full L2
% error and many more.

% The script is divided into multiple steps, based on Bartels adaptivity
% algorithm. Some functions were provided by Kovács, Bartels and Funken

%% Initializing data

TOL = set_TOL; % Error Tol

tau = (T_end-T_0)/16; %First timestepsize

eta_space_global = inf; %First instanciation of the value, is later given by eta_global

%Initialise problem, and store inital conditons into solution space
t = T_0;
j= 1;
f_h = f(Nodes_0,T_0);
u_old = u_0(Nodes_0);

OUTPUT(1,:) = {Nodes_0,Elements_0,Dirichlet_Edges_0,u_old,u_old,T_0,0, 0,f_h,zeros(size(Nodes_0,1)),zeros(size(Nodes_0,1))};

%Nodes_0 store inital grid
Nodes = Nodes_0; Elements = Elements_0;Dirichlet_Edges = Dirichlet_Edges_0;
neumann = [];

while t < T_end
    t = t + tau; %New timestep
    
    while eta_space_global > TOL       
        %% Step 2: Solve for u_new on given grid
        % We will assemble the scheme relevant matrices to determine u_new
        [M,A] = surface_assembly_P1(Nodes,Elements); % Mass and stiffness matrices
        % The implicit scheme gives (M+\tau A) u_j = M (u_(j-1)+ \tau *f_func_j)
        f_h = f(Nodes,t);
        u_new = solve_SFEM(u_old, Nodes, Dirichlet_Edges, M, A, t, tau, f_h);
        
        %% Determine spatial error component
        eta_space_sq = space_estimator(Nodes,Elements,Dirichlet_Edges,u_new,u_old,t, tau, f_h); %Berechnet den Lokalen Fehlerindikator (jeder Eintrag steht für den Fehler auf dem Element)
        eta_space_global = sum(eta_space_sq); %Gesamtfehler

        if eta_space_global < TOL
            break
        end
        
        %% Step 2.A: Mark& Refine
        if marking_criterion == 1
            marked = ( eta_space_sq > theta_refine * max(eta_space_sq) );
        elseif marking_criterion == 2    
            [sorted,pos_mark] =sort(eta_space_sq);
            j = 0; sum_eta = 0;
            while sum_eta < (1-theta_refine)*eta_space_global
                sum_eta = sum_eta+sorted(end-j);
                j = j+1;
            end
            logical_array = zeros(1,size(eta_space_sq,1));
            logical_array(pos_mark((end-j):end)) = 1;
            marked = logical(logical_array);
        end
        
        %% Refine by RGB or NVB
        if refinement_type == 1
            [Nodes,Elements,Interpolation_Nodes, Dirichlet_Edges] = TrefineRGB(Nodes,Elements,Dirichlet_Edges,neumann,marked);
        elseif refinement_type == 2
            [Nodes,Elements,Interpolation_Nodes, Dirichlet_Edges] = refineNVB(Nodes,Elements,Dirichlet_Edges,neumann,marked);
        else
            error(['Refinement Type not correctly initialised'])
        end   
                
        %% Step 2.B: Lift new set of nodes onto Gamma
        
        Nodes = lift(Nodes,d);
        
        %% Step 2.C: Determine nodal values of u_old of the new mesh
        % Evaluation of u_old at new nodes
        if j == 1 %For the first step we know u_old = u_0 exactly, use this information for intialization
            u_old = u_0(Nodes);
        else    
            u_int = 0.5 .* sum(u_old(Interpolation_Nodes(:,1:2)),2); %Take the mean of the two neighbouring nodes, this works due to linear basis functions
            u_old = [u_old; u_int]; %Concatenate the new interpolated nodes to the old function
        end
                       
    end

%% Compute Temporal Error
    eta_temp_sq = temp_estimator(Nodes,Elements,Dirichlet_Edges,u_new, u_old ,t,tau,f);
    eta_temp_global = sum(eta_temp_sq);

%% Step 3: Check if Temporal Error under TOL and store solution
    if eta_temp_global <= TOL
        if j == 1 %Overwrite the initial mesh with the one refined in the first step
            OUTPUT(1,:) = {Nodes,Elements,Dirichlet_Edges,u_old,u_old,T_0,0, 0,f_h,[],[]};
        end
        %% Step 3.A
        tau = 2*tau;
        if t+tau > T_end
            tau = T_end-t;
        end
        j = j+1;
        %Speicher Lösung und Quantities zum Plotten in einem Cell array
        OUTPUT(j,:) = {Nodes,Elements,Dirichlet_Edges,u_new, u_old ,t,eta_space_global, eta_temp_global, f_h, M, A}; %Also store u_old which is already nodally interpolated to be able to write the solutions for continous t values 
        
        %Update das die neue Lösung jetzt die alte ist
        u_old = u_new;
        
        disp(['Part of solution found at t = ', num2str(t)])
        
        %% Step 4: Coarsening       
         while eta_space_global < TOL
            %theta_coarse = 0.25; %Marking strategy parameter
            N0 = size(Nodes_0,1);
            eta_space = space_estimator(Nodes,Elements,Dirichlet_Edges,u_new,u_old,OUTPUT{j,6}, tau/2, f_h); %Berechnet den Lokalen Fehlerindikator (jeder Eintrag steht für den Fehler auf dem Element)
            eta_space_global = sum(eta_space); %Gesamtfehler
            
            %% Marking for Coarsening
            if marking_criterion == 1          
                marked = ( eta_space < theta_coarse * max(eta_space) );
            elseif marking_criterion == 2    
                [sorted,pos_mark] =sort(eta_space);
                k = 1; sum_eta = 0;
                while sum_eta < theta_coarse*eta_space_global
                    sum_eta = sum_eta+sorted(k);
                    k = k+1;
                end
                logical_array = zeros(1,size(eta_space,1));
                logical_array(pos_mark((end-k):end)) = 1;
                marked = logical(logical_array);
            end
            Nodes_before_coarsen = Nodes;

            %% Check which Coarsening should be done
            if coarsen_type == 0 %Naiv
                break
            elseif coarsen_type == 1
                if coarsen_type ~= refinement_type
                   error(['Coarsen Type not the same as Refinement Type'])
                end
                [Nodes,Elements,index_coarse_nodes] = TcoarsenRGB(N0,Nodes,Elements,Dirichlet_Edges,neumann,marked);
            elseif coarsen_type == 2
                if coarsen_type ~= refinement_type
                   error(['Coarsen Type not the same as Refinement Type'])
                end
                [Nodes,Elements,index_coarse_nodes] = coarsenNVB(N0,Nodes,Elements,Dirichlet_Edges,neumann,marked);
            elseif coarsen_type == 3
                Nodes = Nodes_0;
                Elements = Elements_0;
                index_coarse_nodes = 1:size(Nodes_0,1);
            else
                error(['Coarsen Type not correctly initialized'])
            end
            Nodes_after_coarsen = Nodes;
            if size(Nodes_before_coarsen,1) == size(Nodes_after_coarsen,1)
                break %Break if no nodes where coarsend even though the error might not be under the tolerance
            end
            %% Step 4.A
            u_old = u_old(index_coarse_nodes);
            u_new = u_new(index_coarse_nodes);
            f_h = f_h(index_coarse_nodes);
            %% Step 4.B
            eta_space = space_estimator(Nodes,Elements,Dirichlet_Edges,u_new,u_old,OUTPUT{j,6}, tau/2,f_h); %Berechnet den Lokalen Fehlerindikator (jeder Eintrag steht für den Fehler auf dem Element)
            eta_space_global = sum(eta_space); %Gesamtfehler
         end
        eta_space_global = inf; %Reset Error quantaties to get back in the while loops
    else
        %% Step 3.B
        t = t-tau; %Set the time back to u_old to Resolve u_new with half timestepsize to get smaller eta_space error.
        tau = tau/2; %New (halfed) timestepsize
        eta_space_global = inf; %Reset Error quantaties to get back in the while loops
    end       
end
end