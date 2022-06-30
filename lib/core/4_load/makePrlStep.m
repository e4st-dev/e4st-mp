function gencost = makePrlStep(mpc, esc, hour)
%makePrlStep: Create gencost steps for price-responsive load in base hour
%
% E4ST
% Copyright (c) 2012-2018 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Get the parameters
elsty = esc.elasticity;
dist_cost = esc.dist_cost * (1 - esc.dist_loss); % $ / Wholesale MWh
price_vec = esc.prl_vec;
max_price_vec = esc.max_prl_vec;

% distribution cost is per load + losses
% uplift is retail MWh

%% Calculate the pivot points and default prices
%bus_table = array2table([mpc.bus(:, 1), res_bus.load_hourly(:, hour), res_bus.lmp_hourly(:, hour)], 'VariableNames', {'bus', 'load', 'lmp'});
bus_table = array2table([mpc.bus(:, 1), mpc.default_load(:, hour), mpc.default_price(:, hour)], 'VariableNames', {'bus', 'load', 'lmp'});
%bus_table = array2table([mpc.bus(:, 1), mpc.hourly_load(:, hour), mpc.default_price(:, hour)], 'VariableNames', {'bus', 'load', 'lmp'});
gen_table = join(array2table(mpc.gen(:, 1), 'VariableNames', {'bus'}), bus_table);
default_load = gen_table{:, 'load'}; % MW
default_price = gen_table{:, 'lmp'}; % $ / Wholesale MWh

pivot_load = default_load; % MW
pivot_price = default_price + dist_cost; % Add the distribution cost to pivot price

%% Apply load growth to pivot loads (from the first simulation year)
info = esc.load_growth;
info = filterInfo(info, esc);

info.value = info.value(:, ['Y', num2str(esc.year)]);

for i = 1:height(info.value)
    idx_info = i == (1:height(info.value))';
    cur_info = filterStruct(info, idx_info);
    idx_gen = getInfoIdx(mpc, esc, cur_info, 'gen');

    if ~any(idx_gen)
        continue
    end

    load_growth = cur_info.value{:, :}.^(esc.year - esc.first_sim_yr);
    pivot_load(idx_gen) = pivot_load(idx_gen) * load_growth;
end

%% Select 'dl' only - Find loads needed to add step gencost
idx_dl = strcmp(mpc.genfuel(:, 1), 'dl');
n_dl = sum(idx_dl);
pivot_load = pivot_load(idx_dl);
pivot_price = pivot_price(idx_dl);

if isfield(mpc, 'uplift')
    uplift = mpc.uplift(idx_dl); % Should be $ / Wholesale MWh
else
    uplift = zeros(n_dl, 1);
end

%% Create PRL Steps
% Transform input vector of price percentages into actual prices
% For the last two steps, set very high prices, unless the percentages
% are actually greater than these high prices ($5,000 and $10,000)
% Use 10 steps of prices
step_price = repmat(pivot_price, 1, 10) .* price_vec';
step_price(:, 1:8) = min(step_price(:, 1:8), max_price_vec(1:8)'); % ceiling
step_price(:, 9:10) = max(step_price(:, 9:10), max_price_vec(9:10)'); % floor

% For the vertexPower, calculate the amount of power demanded
% at each price in the step_price
% Set the power demanded at the top price level equal to 0
% Calculate the actual vertex price and power
vertex_power = zeros(n_dl, 10);
vertex_price = step_price;
for i = 1:9
    vertex_price(:, i) = mean([step_price(:, i), step_price(:, i + 1)], 2);
    vertex_power(:, i) = pivot_load .* ((step_price(:, i)) ./ (pivot_price)).^elsty;
    if ~isreal(vertex_power(:, i))
        fprintf('Huge negative price at hour %d\n', hour);
    end
end
vertex_power(:, 10) = 0;
vertex_power = -vertex_power; % 'dl' have negative sign

% Calculate the total costs, or f(p) in the gencost notation
% Need to do this iteratively, as each (Power*Price) block must be added
% to the result of the previous TotalCostVector block
total_cost_vec = zeros(n_dl, 10);
for i = 9:-1:1
    % Subtract distribution cost and uplift (if applicable) from vertexPrice
    total_cost_vec(:, i) = total_cost_vec(:, i+1) + ...
        (vertex_power(:, i) - vertex_power(:, i + 1)) .* (vertex_price(:, i) - dist_cost - uplift);
end

%% Export the resultant gencost: p(MW), f(p) ($/hr)
gencost = zeros(n_dl, 24);
gencost(:, 1:4) = repmat([1, 0, 0, 10], n_dl, 1);
gencost(:, 5:2:23) = vertex_power(:, 1:10);
gencost(:, 6:2:24) = total_cost_vec(:, 1:10);

% Filter out and reset gencost with negative LMP and zero loads, which result in complex power
idx_reset = (default_price(idx_dl) <= 0) | (default_price(idx_dl) >= 5000) | (default_load(idx_dl) <= 0);
if any(idx_reset)
    gencost(idx_reset, 1:5) = repmat([2, 0, 0, 2, 5000], sum(idx_reset), 1);
    gencost(idx_reset, 6:end) = 0;
end

% Check gencost
if ~isreal(gencost)
    fprintf('Warning: Unrealistic gencosts at hour %d\n', hour);
    pause();
end
