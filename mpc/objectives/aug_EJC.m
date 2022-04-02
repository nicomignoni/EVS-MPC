function [sol, J] = aug_EJC(k, H, P, c, p, Max, rho, l, d_opt, x)
    prob = optimproblem("ObjectiveSense", "min");
    
    r = optimvar("r", 1, H, "LowerBound", 0, "UpperBound", Max.r);
    d = optimvar("d", P, H, "LowerBound", 0, "UpperBound", Max.d);

    x0.d = squeeze(x.P.d(:,k,:));

    prob.Constraints.OptChar = r + sum(d) == d_opt(k,:);      
    
    prob.Objective = EJC(k, d, r, c, p, l) ...
                   + sum(0.5*rho.C*(d - x0.d).^2, "all");
                  
    options = optimoptions("quadprog", "Display", "off"); 
    sol = solve(prob, "options", options);
    J = evaluate(EJC(k, d, r, c, p, l), sol);
end

