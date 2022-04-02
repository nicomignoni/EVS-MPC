function [sol, J] = aug_EJP(i, H, C, M, G, D, c, p, Max, alpha, eta, rho, init, l, x)
    prob = optimproblem("ObjectiveSense", "min");
    
    p_up   = optimvar("p_up", 1, H, "LowerBound", 0, "UpperBound", Max.p);
    p_down = optimvar("p_down", 1, H, "LowerBound", 0, "UpperBound", Max.p);
    q_up   = optimvar("q_up", M, H, "LowerBound", 0, "UpperBound", Max.q);
    q_down = optimvar("q_down", M, H, "LowerBound", 0, "UpperBound", Max.q);
    d      = optimvar("d", C, H, "LowerBound", 0, "UpperBound", Max.d);
    e      = optimvar("e", M, H, "LowerBound", 0);
    
    x0.q_up   = squeeze(x.M.q_up(i,:,:));   
    x0.q_down = squeeze(x.M.q_down(i,:,:));
    x0.d      = squeeze(x.C.d(i,:,:)); 
    x0.e      = squeeze(x.M.e(i,:,:));
    
    prob.Constraints.Balance = G(i,:) - D(i,:) + p_up - p_down - sum(d) ...
                             + sum(q_down - q_up) == zeros(1,H);
    prob.Constraints.Storage = e == alpha*p.M.*[init.e(i,:)' e(:,1:end-1)] ...
                             + eta.chr.*q_up - eta.dis.*q_down;
                              
    prob.Objective = EJP(i, p_up, p_down, q_up, q_down, d, e, p, c, l) ...
                   + sum(0.5*rho.M*(q_up - x0.q_up).^2, "all") ...
                   + sum(0.5*rho.M*(q_down - x0.q_down).^2, "all") ...
                   + sum(0.5*rho.C*(x0.d - d).^2, "all") ...
                   + sum(0.5*rho.M*(x0.e - e).^2, "all");
                   
    options = optimoptions("quadprog", "Display", "off"); 
    sol = solve(prob, "options", options);
    J  = evaluate(EJP(i, p_up, p_down, q_up, q_down, d, e, p, c, l), sol); 
end

