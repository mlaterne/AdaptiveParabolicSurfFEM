function nodes = lift(x,d,grad_h)
%This function takes a set of nodes x on Gamma_h and approximates a(x), its
%corresponding point on Gamma via surface normal projection.

% Elements and Dirichlet_Edges will stay the same only the Nodes change
% coordiantes!

% If the distance function is known we use the distmesh implementation for
% the lift (i.e. computing the gradient with finite difference) (grad_d is
% not required!)

% If we only know a sufficiently smooth level set function h:
% Use a first order approximation of the lifts as described in Thesis

switch nargin
    case 3
    h = d;
    nodes = zeros(size(x,1),3);
    TOL = 5*10^-10;
    err = 1;
    %l = logical(1 == ones(size(x,1),1)); sollte mit logical recht einfach
    %vektorisierbar sein
    for i = 1:size(x,1)
        a = h(x(i,:)); b = grad_h(x(i,:)); %Inital grad and distance function
        if abs(a) < TOL
            nodes(i,:) = x(i,:);
            continue
        end
        x_new = x(i,:); %Intermediate x
        while err > TOL
            x_tilde = x_new-(a*b)./(norm(b)^2); %One approximation step
            dist = sign(h(x(i,:)))*norm(x_tilde- x(i,:));
            x_new = x(i,:)- dist.* grad_h(x_tilde)./norm(grad_h(x_tilde));
            a = h(x_new); b = grad_h(x_new); %Determine new distance and gradient
            % Due to probably parametrization differences i need to take the
            % negative of -x_new+x(i,:) to obtain correct results
            err = (a^2/norm(b)^2 + norm(b/norm(b) +sign(h(x(i,:)))* (x_new-x(i,:))./norm(x_new-x(i,:)))^2)^1/2; %Second term ensure that the normal direction from the surface points to x_0
            %disp(['Distance function = ', num2str(a), ' and error = ', num2str(err)])
        end
        nodes(i,:) = x_new;
        err = 1;
    end
    case 2
        h0 = 0.1;
        deps=sqrt(eps)*h0;
        TOL = 5*10^-10;
        err = 1;
        d_nodes=d(x);
        while err > TOL
            dgradx=(d(x+deps.*ones(1,size(x,1))'*[1 0 0])-d_nodes)/deps; % Numerical
            dgrady=(d(x+deps.*ones(1,size(x,1))'*[0 1 0])-d_nodes)/deps; % gradient
            dgradz=(d(x+deps.*ones(1,size(x,1))'*[0 0 1])-d_nodes)/deps; %
            dgrad2=(dgradx.^2+dgrady.^2+dgradz.^2);
            nodes = x-[d_nodes.*dgradx,d_nodes.*dgrady,d_nodes.*dgradz]./dgrad2;
            x = nodes;
            d_nodes = d(nodes);
            err = max(d_nodes);
            
        end
end
end
