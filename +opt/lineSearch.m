function out = lineSearch(x0, s, f, t0, steps, stepMod)
    % SETUP
    fCalls = 0;
    z = f(x0); fCalls = fCalls+1;
    x = zeros(length(x0), 1); 
    x(:, end) = x0;           % Enforcing column vectors
    t = t0;
    iter = 1; iterTot = 1;
    
    % ITERATING UNTIL MINIMUM BRACKETED
    while length(z) < 3 || z(end) < z(end-1)
        % INCREASING STEP SIZE IF NOT BRACKETING
        if iter == steps
            t = stepMod*t;
            iter = 0;
        end
        
        % STEPPING AND EVALUATING
        x = cat(2, x, x(:, end) + t*s);
        z(end+1) = f(x(:, end)); fCalls = fCalls+1;
        iter = iter+1;
        iterTot = iterTot+1;
        
        if iter == 31
            pts = [];
            warning("lineSearch() ran out of iterations")
            return
        end
    end
    
    % OUTPUTTING POINTS
    pts = struct('x', x(:, end-2), 'z', z(end-2));
    pts(end+1) = struct('x', x(:, end-1), 'z', z(end-1));
    pts(end+1) = struct('x', x(:, end), 'z', z(end));
    out.pts = pts;
    out.iters = iterTot;
    out.allPts = x;
    out.fCalls = [fCalls 0 0];
end