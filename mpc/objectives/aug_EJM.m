function [sol, J] = aug_EJM(j, H, P, c, p, Max, alpha, eta, zeta, rho, init, l, x)
    prob = optimproblem("ObjectiveSense", "min");
    
    q_up   = optimvar("q_up", P, H, "LowerBound", 0, "UpperBound", Max.q);
    q_down = optimvar("q_down", P, H, "LowerBound", 0, "UpperBound", Max.q);
    e      = optimvar("e", P, H, "LowerBound", 0);
    
    x0.q_up   = squeeze(x.P.q_up(:,j,:));
    x0.q_down = squeeze(x.P.q_down(:,j,:));
    x0.e      = squeeze(x.P.e(:,j,:));
    
    prob.Constraints.BatteryLimit = sum(e) <= Max.b.M(j) - init.b.M(j);
    prob.Constraints.Battery = e == alpha*ones(P,1).*p.M(j,:).*[init.e(:,j) e(:,1:end-1)] ...
                                  + eta.chr.*q_up - eta.dis.*q_down;
                              
    prob.Objective = EJM(j, P, q_up, q_down, e, p, zeta, c, l) ...
                   + sum(0.5*rho.M*(x0.q_up - q_up).^2, "all") ...
                   + sum(0.5*rho.M*(x0.q_down - q_down).^2, "all") ...
                   + sum(0.5*rho.M*(e - x0.e).^2, "all");
              
    options = optimoptions("quadprog", "Display", "off"); 
    sol = solve(prob, "options", options);
    J   = evaluate(EJM(j, P, q_up, q_down, e, p, zeta, c, l), sol);
    
%     % Solved the J-constraint version in case the there is an economic loss
%     if J > 0
% %         xi.q_up = squeeze(x0.M.q_up(:,j,:));
% %         xi.q_down = squeeze(x0.M.q_down(:,j,:));
% %         xi.e = squeeze(x0.M.e(:,j,:));
%     
%         options = optimoptions("fmincon", "Display", "off"); 
%         prob.Constraints.Revenue = EJM(j, P, q_up, q_down, e, p, zeta, c, l) <= 0;
%         sol = solve(prob, sol, "options", options);
%         J   = evaluate(EJM(j, P, q_up, q_down, e, p, zeta, c, l), sol);
%     end
end

