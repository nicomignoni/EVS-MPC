function [M,Max,init,zeta,p] = update_M(t,T,H,Max,init,ev_battery,zeta,p)
    n_sample = size(ev_battery.Capacity,1);
    M = size(p.M,1);
    P = size(init.e,1);
    still_parked = rand(M,1) <= p.M(:,t);

    % Remove max, and init capacity and storage
    Max.b.M  = Max.b.M(still_parked);
    init.b.M = init.b.M(still_parked);
    init.e = init.e(:,still_parked);

    % Remove probability
    p.M = p.M(still_parked,:);

    % Remove degradation coefficient 
    zeta = zeta(still_parked);

    for j=1:poissrnd(1)
        % Create a random index
        index = randsample(n_sample, 1);
        
        % Create probability 
        mu    = 3*randi([1, 24], 1, 1);
        sigma = randi([1, 24], 1, 1);
        pd  = makedist("Normal", "mu", mu, "sigma", sigma);
        tpd = truncate(pd, 0, inf); 
        p.M = [p.M; zeros(1,t-1) 1-cdf(tpd, 1:24*T+H-t+1)];

        % Add new max capacity, intial capacity and degradation coefficient 
        Max.b.M = [Max.b.M; ev_battery.Capacity(index)];
        init.b.M = [init.b.M; 0];
        init.e = [init.e zeros(P,1)];
        zeta = [zeta 1e-2*rand()];
    end
    M = size(p.M,1);
end

