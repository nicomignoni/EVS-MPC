% General Settings
clear
beep off;
rng('default');

%% Figure settings
% Interpreter settings
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultTextarrowshapeInterpreter','latex');

% Figure settings
set(groot, 'DefaultFigureUnits', "centimeters");
set(groot, 'DefaultFigurePaperType', '<custom>');
set(groot, 'DefaultFigurePaperSize', [9 6]);
set(groot, 'DefaultFigurePosition', [12 4 9 6]);
set(groot, 'DefaultFigurePaperPositionMode', 'auto');

% Axes fontsize
set(groot, 'DefaultAxesFontSize', 8);
% set(groot, 'DefaultAxesUnits', 'normalized');
% set(groot, 'DefaultAxesOuterPosition', [0 0 1 0.95]);
% set(groot, 'DefaultAxesInnerPosition', [0.13 0.18 0.72 0.7]);

%% Loading data
load("data/energy_data.mat");
load("data/ev/ev_battery.mat")

%% Constant settings
T = 3;     % Simulation time (days)
H = 24;     % Control horizon (h)
C = 10;    % Number of charging EVs
M = 10;    % Number of trading EVs
P = 5;    % Number of prosumers

% Energy flow upper bounds
Max.p = 10;
Max.q = 10;
Max.d = 10;
Max.r = 100;

% Prices from Enel Bioraria fees table
% (https://www.enel.it/it/luce-e-gas/luce/offerte/e-light-bioraria)
price.blue   = 0.13599;
price.orange = 0.15803;
epsilon = 0.7;                    % Scale factor between C_sell and C_buy
c.sell = [price.blue*ones(1,8) ...
          price.orange*ones(1,11) ... 
          price.blue*ones(1,5)];  % Retailer -> Users price
c.buy  = epsilon*c.sell;          % Users -> Retailer price

% ESSs leakage and efficiencies
eta.chr = 0.98;
eta.dis = 1.02;
alpha   = 1 - 1e-2;

% Prosumers' demand and generation
n_sample = size(DP, 1);
D = [repmat(DP, [1 T]) DP(:,1:H)];
G = [repmat(PV, [1 T]) DP(:,1:H)];
D = D(randsample(n_sample, P),:) + 3*rand(P,24*T + H);
G = G(randsample(n_sample, P),:) + 3*rand(P,24*T + H);
new_max_D = 3;  
max_D = max(D,[],"all");
min_D = min(D,[],"all");
D = min_D + (new_max_D - min_D).*(D - min_D)./(max_D - min_D);

%% Variable settings
% tEvs storage parameters
zeta = 1e-2*rand(1,M);

% Capacity Data 
% (https://ev-database.org/cheatsheet/range-electric-car)
ni      = 1 - 1e-2;
fc_frac = 0.3;
index.C = randi(size(ev_battery.Capacity,1), 1, C);
index.M = randi(size(ev_battery.Capacity,1), 1, M);

Max.b.C  = ev_battery.Capacity(index.C); 
Max.b.M  = ev_battery.Capacity(index.M);
init.b.C = 0*rand(C,1).*Max.b.C; 
init.b.M = 0.5*rand(M,1).*Max.b.M; 
init.e   = zeros(P,M);

h_rc  = ev_battery.Range(index.C)./ev_battery.ChargeSpeed(index.C);
gamma = -h_rc.^(-1).*log((ni - 1).*Max.b.C./(init.b.C - Max.b.C));

% Exponential charging profile
b_opt = Max.b.C + (init.b.C - Max.b.C).*exp(-(1:24*T+H).*gamma);
d_opt = min((b_opt - alpha*[init.b.C b_opt(:,1:end-1)])/eta.chr, Max.d + Max.r);
 
% Probabilities
mu.C    = randi([1, 24], 1, C);
mu.M    = 3*randi([1, 24], 1, M);
sigma.C = randi([1, 24], 1, C);
sigma.M = randi([1, 24], 1, M);

for k=1:C
    pd  = makedist("Normal", "mu", mu.C(k), "sigma", sigma.C(k));
    tpd = truncate(pd, 0, inf); 
    p.C(k,:) = 1 - cdf(tpd, 1:24*T+H);
end
for j=1:M
    pd  = makedist("Normal", "mu", mu.M(j), "sigma", sigma.M(j));
    tpd = truncate(pd, 0, inf); 
    p.M(j,:) = 1 - cdf(tpd, 1:24*T+H);
end