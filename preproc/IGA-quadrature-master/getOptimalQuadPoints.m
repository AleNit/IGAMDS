
function [w, x, rec, it] = getOptimalQuadPoints(knot, p, w0, x0, knot0, rec)

warn_state  = warning();
warning off;
n  = numel(knot)-p-1; % dimension of our spline space

% compute all greville points and integrals (used for initial guess)
greville       = zeros(n,1);
exact_integral = zeros(n,1);
for i=1:n
    greville(i)       = sum(knot(i+1:i+p)) / p;
    exact_integral(i) = (knot(i+p+1)-knot(i))/(p+1);
end

% exact_integraln=accurate_numeric(knot,p);
if exist('x0') % if initial guess is provided, use these
    x = x0;
    w = w0;
else           % else compute them based on greville points and integrals
    w   = (exact_integral(1:2:end) + exact_integral(2:2:end))  ;
    x   = (      greville(1:2:end) +       greville(2:2:end))/2;
    rec = 1;     % counter variable to count the number of recursive calls
end

dF = sparse(n,n);
newton_tol       = 1e-15;    % convergence tolerance
newton_max_it    = 150;      % max iterations before divergence
while true                   % recursive loop from algorithm 2
    for it = 1:newton_max_it   % newton iteration loop
        [N,dN] = BSpline(knot, p, x);
        F      = N*w - exact_integral;
        %      dF     = [N, dN*diag(sparse(w))];
        
        dF(:,1:2:n) = N;
        dF(:,2:2:n) = dN*diag(sparse(w));
        
        dx = dF \ -F;
        w = w + dx(1:2:end);
        x = x + dx(2:2:end);
        
        % test for diverging (coarse heuristic, see section 3.3)
        if( min(x)<knot(1)  )     
            break;   
        end
        if( max(x)>knot(end))     
            break;   
        end
        
        % test for converging
        if(norm(dx)<newton_tol)  
            warning(warn_state); 
            return;  
        end
    end
    
    % at this point, newton iteration has diverged. solve recursively on easier knot
    if exist('knot0')
        [w, x, rec] = getOptimalQuadPoints((knot0 + knot)/2, p, w0, x0, knot0, rec);
        knot0       = (knot0 + knot)/2;
    else
        uniformKnot = linspace(knot(1),knot(end), n+p+1);
        [w, x, rec] = getOptimalQuadPoints(uniformKnot, p);
        knot0       = uniformKnot;
    end
    rec = rec + 1;
    x0 = x;
    w0 = w;
end % loop up and start newton iteration with better initial guess

end


function exact_integral=accurate_numeric(knot,p)

ng=8;
[weightv,coordv]=quadrature(ng, 'GAUSS', 1 );
uknot=unique(knot);
ig=1;
for i=2:length(uknot)
    xiE=[uknot(i-1),uknot(i)];
    J2      = 0.5*(xiE(2)-xiE(1));
    Xi(ig:ig+ng-1) = 0.5 * ( ( xiE(2) -xiE(1) ) * coordv(:,1) + xiE(2) + xiE(1));
    W(ig:ig+ng-1)=weightv*J2;
    ig=ig+ng;
end
[N ~] = BSpline(knot, p, Xi);
exact_integral= N*W';

end
