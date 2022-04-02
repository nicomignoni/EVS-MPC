function [C,Max,init,d_opt,p] = update_C(t,T,H,Max,init,ev_battery,ni,alpha,eta,d_opt,p)
    n_sample = size(ev_battery.Capacity,1);
    C = size(p.C,1);
    still_parked = rand(C,1) <= p.C(:,t);
    
    % Remove max and init capacity
    Max.b.C  = Max.b.C(still_parked);
    init.b.C = init.b.C(still_parked);

    % Remove probability
    p.C = p.C(still_parked,:);

    % Remove optimal intake
    d_opt = d_opt(still_parked,:);

    for k=1:poissrnd(1)
        % Create a random index
        index = randsample(n_sample, 1);

        % Add new max and intial capacity
        Max.b.C = [Max.b.C; ev_battery.Capacity(index)];
        init.b.C = [init.b.C; 0];

        % Add probability distribution
        mu    = randi([1, 24], 1);
        sigma = randi([1, 24], 1);
        pd  = makedist("Normal", "mu", mu, "sigma", sigma);
        tpd = truncate(pd, 0, inf);
        p.C = [p.C; zeros(1,t-1) 1-cdf(tpd, 1:24*T+H-t+1)];

        % Add recharge time and gamma
        h_rc = ev_battery.Range(index)./ev_battery.ChargeSpeed(index);
        gamma = -1/h_rc.*log((ni - 1).*Max.b.C(end)./(init.b.C(end) - Max.b.C(end)));

        % Add new optimal intake curve
        b_opt = Max.b.C(end) + (init.b.C(end) - Max.b.C(end)).*exp(-(1:24*T+H-t+1).*gamma);
        d_opt = [d_opt; 
                zeros(1,t-1) min((b_opt - alpha*[init.b.C(end) b_opt(:,1:end-1)])/eta.chr, Max.d + Max.r)];
    end
    C = size(p.C,1);
end

