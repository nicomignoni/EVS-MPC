function [x,residual,l] = ADALM(P,C,M,H,p,c,d_opt,G,D,Max,init,alpha,eta,zeta,ITER,chi,rho)
    %%%%%%% Accelerated Distributed Agumented Lagrangian Method %%%%%%%
    addpath("./objectives");
    % Intial solution (iteration 0)
    % (Charging Evs)
    x.C.r   = zeros(C,H);
    x.C.d   = zeros(P,C,H);
    x.C.b   = zeros(C,H);
    x.C.J   = zeros(1,C); 
    % (Trading Evs)
    x.M.q_up   = zeros(P,M,H);
    x.M.q_down = zeros(P,M,H);
    x.M.e      = zeros(P,M,H);
    x.M.J      = zeros(1,M); 
    % (Prosumers)
    x.P.p_up   = zeros(P,H);
    x.P.p_down = zeros(P,H);
    x.P.q_up   = zeros(P,M,H);
    x.P.q_down = zeros(P,M,H);
    x.P.d      = zeros(P,C,H);
    x.P.e      = zeros(P,M,H);
    x.P.J      = zeros(1,P); 

    % Lambda
    l.psi      = zeros(P,C,H);
    l.phi_up   = zeros(P,M,H);
    l.phi_down = zeros(P,M,H);
    l.theta    = zeros(P,M,H);
    
    % Residuals
    residual.q_up.all   = zeros(P,M,H);
    residual.q_down.all = zeros(P,M,H);
    residual.d.all      = zeros(P,C,H);
    residual.e.all      = zeros(P,M,H);
    
    vn = fieldnames(residual);
    for i=1:numel(vn)
        residual.(vn{i}).mean = zeros(1,ITER);
        residual.(vn{i}).std  = zeros(1,ITER);
    end
    
    % Distributed ALM
    for tau=1:ITER
       tic;
       for k=1:C
           [sol, J] = aug_EJC(k, H, P, c, p, Max, rho, l, d_opt, x);
            
           x.C.r(k,:)   = x.C.r(k,:) + chi*(sol.r - x.C.r(k,:));
           x.C.d(:,k,:) = squeeze(x.C.d(:,k,:)) + chi*(sol.d - squeeze(x.C.d(:,k,:)));
           x.C.J(k)     = J;
       end
       
       for j=1:M
           [sol, J] = aug_EJM(j, H, P, c, p, Max, alpha, eta, zeta, rho, init, l, x);
           
           x.M.q_up(:,j,:)   = squeeze(x.M.q_up(:,j,:)) + chi*(sol.q_up - squeeze(x.M.q_up(:,j,:)));
           x.M.q_down(:,j,:) = squeeze(x.M.q_down(:,j,:)) + chi*(sol.q_down - squeeze(x.M.q_down(:,j,:)));
           x.M.e(:,j,:)      = squeeze(x.M.e(:,j,:)) + chi*(sol.e - squeeze(x.M.e(:,j,:)));
           x.M.J(j)          = J;
       end
       
       for i=1:P
           [sol, J] = aug_EJP(i, H, C, M, G, D, c, p, Max, alpha, eta, rho, init, l, x);
           
           x.P.p_up(i,:,:)   = x.P.p_up(i,:,:) + chi*(sol.p_up - x.P.p_up(i,:,:));
           x.P.p_down(i,:,:) = x.P.p_down(i,:,:) + chi*(sol.p_down - x.P.p_down(i,:,:));
           x.P.q_up(i,:,:)   = squeeze(x.P.q_up(i,:,:)) + chi*(sol.q_up - squeeze(x.P.q_up(i,:,:)));
           x.P.q_down(i,:,:) = squeeze(x.P.q_down(i,:,:)) + chi*(sol.q_down - squeeze(x.P.q_down(i,:,:)));
           x.P.d(i,:,:)      = squeeze(x.P.d(i,:,:)) + chi*(sol.d - squeeze(x.P.d(i,:,:)));
           x.P.e(i,:,:)      = squeeze(x.P.e(i,:,:)) + chi*(sol.e - squeeze(x.P.e(i,:,:)));
           x.P.J(i)          = J;
       end
       
       % Update residuals
       residual.q_up.all   = x.P.q_up - x.M.q_up; 
       residual.q_down.all = x.P.q_down - x.M.q_down;
       residual.d.all      = x.C.d - x.P.d;
       residual.e.all      = x.M.e - x.P.e;
       
       % Update multipliers
       l.psi      = l.psi + chi*rho.C.*residual.d.all;
       l.theta    = l.theta + chi*rho.M.*residual.e.all;
       l.phi_up   = l.phi_up + chi*rho.M.*residual.q_up.all;
       l.phi_down = l.phi_down + chi*rho.M.*residual.q_down.all;
       
       % Add residuals mean and std 
       residual.q_up.mean(tau)   = mean(residual.q_up.all, "all");
       residual.q_up.std(tau)    = std(residual.q_up.all, 0, "all");
       residual.q_down.mean(tau) = mean(residual.q_down.all, "all");
       residual.q_down.std(tau)  = std(residual.q_down.all, 0, "all");
       residual.d.mean(tau)      = mean(residual.d.all, "all");
       residual.d.std(tau)       = std(residual.d.all, 0, "all");
       residual.e.mean(tau)      = mean(residual.e.all, "all");
       residual.e.std(tau)       = std(residual.e.all, 0, "all");
       
       time_elapsed = toc;
       
       % Iteration summary
       fprintf("Iter. %d completed. Elapsed time: %.2f (s) \n", tau, time_elapsed);
       fprintf("Res. q_up: %d\n", residual.q_up.mean(tau));
       fprintf("Res. q_down: %d\n", residual.q_down.mean(tau));
       fprintf("Res. d: %d\n", residual.d.mean(tau));
       fprintf("Res. e: %d\n\n", residual.e.mean(tau));
    end
end