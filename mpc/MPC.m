%% Load settings
clear
beep off
addpath("../utils");
addpath("./agents");
addpath("./objectives");
run("../general_settings.m");
approach = 'mpc';

%% Convergence settings
rho.M = 0.5;
rho.C = 10;
chi   = 0.4;
ITER  = 100;

%% MPC (ADALM)
for t=1:24*T
    % Take prosumer's generation and demand for the control horizon
    Gh = G(:,t:t+H-1);
    Dh = D(:,t:t+H-1);

    % Take probabilities for the control horizon
    ph.C = p.C(:,t:t+H-1);
    ph.M = p.M(:,t:t+H-1);

    % Take optimal energy intake
    dh_opt = d_opt(:,t:t+H-1);
        
    % Solve ADALM
    fprintf("Hour n. %d, C=%d, M=%d, starting ADALM...\n", t,C,M);
    [x_t,residual_t,l_t] = ADALM(P,C,M,H,ph,c,dh_opt,Gh,Dh,Max,init,alpha,eta,zeta,ITER,chi,rho);
    x(t) = x_t; 
    residual(t) = residual_t;
    l(t) = l_t;
    fprintf("ADALM finished.\n");

    % Update initial storage
    init.e = x(t).M.e(:,:,1);

    % Update retailer's energy prices
    c.sell = circshift(c.sell, -1);
    c.buy  = circshift(c.buy, -1);

    % Update cEVs
    [C,Max,init,d_opt,p] = update_C(t,T,H,Max,init,ev_battery,ni,alpha,eta,d_opt,p);
    
    % Update tEVs
    [M,Max,init,zeta,p] = update_M(t,T,H,Max,init,ev_battery,zeta,p);
end

%% Collect data
% Variables and projections on H
fn = [{'p_up'} {'p_down'} {'q_up'} {'q_down'} {'d'} {'e'}];
for i=1:numel(fn)
    projections.(fn{i}) = zeros(H,24*T);
    for t=1:24*T
        for h=1:H
            if ismatrix(x(t).P.(fn{i}))
                projections.(fn{i})(h,t+h-1) = mean(x(t).P.(fn{i})(:,h));
            else 
                projections.(fn{i})(h,t+h-1) = mean(sum(x(t).P.(fn{i})(:,:,h),2),1);
            end
        end
    end
    projections.(fn{i}) = projections.(fn{i})(:,1:24*T);
end

% Energy composition
pros_energy.in  = max(mean(G(1:24*T)) + projections.p_up(1,:) + projections.q_down(1,:), 0);
pros_energy.out = max(mean(D(1:24*T)) + projections.p_down(1,:) + projections.q_up(1,:) + projections.d(1,:), 0);

comp.p_down = projections.p_down(1,:)./pros_energy.out;
comp.d      = projections.d(1,:)./pros_energy.out;
comp.q_up   = projections.q_up(1,:)./pros_energy.out;
comp.D      = mean(D(1:24*T))./pros_energy.out;

comp.p_up   = projections.p_up(1,:)./pros_energy.in;
comp.q_down = projections.q_down(1,:)./pros_energy.in;
comp.G      = mean(G(1:24*T))./pros_energy.in;

%% Visualization
% MPC results and projections
fig = figure("Name", "MPC");
tl = tiledlayout(3,1);

fig.Units = "centimeters";
fig.Position = [0 0 24 9];
fig.PaperUnits = "centimeters";
fig.PaperPosition = [0 0 24 9];
fig.PaperSize = [24 12];

nexttile;
plot_projections(projections.p_down, H, T, [0.0431 0.7098 0.5529], [1 1 1], ...
                "Avg. $\check{a}_{ih}$ \ [kWh]");
nexttile;
plot_projections(projections.e, H, T, [0.3294 0.3882 1], [1 1 1], ...
                "Avg. $u_{ijh}$ \ [kWh]");
nexttile;
plot_projections(projections.d, H, T, [1 0.0941 0.0941], [1 1 1], ...
                 "Avg. $p_{ikh}$ \ [kWh]");
xlabel("$t$ [hours]");
saveas(tl, "..\plot\pdf\mpc_plot.pdf");

% Energy Composition 
fig = figure("Name", "Energy Composition");
        fig.Position = [12 4 18 3.5];
        fig.PaperSize = [18 3.5];
        tl = tiledlayout(1,2, "TileSpacing", "compact");
        tl.Padding = "compact";
        nexttile;
        ax = gca;
        b = bar([comp.p_down;  comp.d; comp.q_up; comp.D]', "stacked", ...
                "EdgeAlpha", 0);
        grid on;
        ax.XLabel.String = "$t$ [hours]";
        ax.YLabel.String = "Out-flow fractions (\%)";
        b(1).FaceColor = col.ret.rgb;
        b(2).FaceColor = col.ev.ch.rgb;
        b(3).FaceColor = col.ev.tr.rgb;
        b(4).FaceColor = col.pros.rgb;
        leg = legend('$\check{a}_{ijh}$', '$p_{ikh}$', '$\hat{d}_{ijh}$', '$D_h$',...
               'Orientation', 'horizontal',...
               'Location', 'north');
        uistack(leg, "top")
        
        nexttile;
        ax = gca;
        b = bar([comp.p_up; comp.q_down; comp.G]', "stacked", ...
                "EdgeAlpha", 0);
        grid on;
        ax.XLabel.String = "$t$ [hours]";
        ax.YLabel.String = "In-flow fractions (\%)";
        b(1).FaceColor = col.ret.rgb;
        b(2).FaceColor = col.ev.tr.rgb;
        b(3).FaceColor = col.pros.rgb;
        legend('$\hat{a}_{ijh}$', '$\check{d}_{ijh}$', '$G_h$', ...
               'Orientation', 'horizontal', ...
               'Location', 'north');
        
        saveas(ax, ['../plot/pdf/' filename '/composition.pdf']);
        savefig(fig, ['../plot/fig/' filename '/composition.fig']);

%% Save data
filename = [approach '_' datestr(now,'mm-dd-yyyy_HH-MM-SS')];
save(['../data/results/' filename '.mat']);