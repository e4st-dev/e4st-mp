function res_bus = result_bus(res, mpc, esc, eopt)
%result_bus Summarize E4ST simulation results by bus

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Hourly Results
define_constants;
hours = esc.hrs_map{:, 'hours'};
probability = esc.hrs_map{:, 'probability'};
nh = length(hours);

load_hourly = zeros(size(res.base.bus, 1), nh);
lmp_hourly = zeros(size(res.base.bus, 1), nh);
for i = 1:nh
    if i == 1
        bus = res.base.bus;
        gen = res.base.gen;
    else
        bus = res.cont(i - 1).bus;
        gen = res.cont(i - 1).gen;
    end
    load_hourly(:, i) = total_load(bus, gen, (1:size(bus, 1))');
    lmp_hourly(:, i) = bus(:, LAM_P) / probability(i);
end

res_bus.load_hourly = load_hourly;
res_bus.lmp_hourly = lmp_hourly;
res_bus.hourly_shape = mpc.hourly_shape;
res_bus.uncurtailed_load_hourly = mpc.hourly_load;

%% Annual Results
res_bus.annual = table(mpc.bus(:, 1), 'VariableNames', {'bus'});
res_bus.annual{:, 'load'} = load_hourly * hours;
res_bus.annual{:, 'uncurtailed_load'} = res_bus.uncurtailed_load_hourly * hours;
res_bus.annual{:, 'curtailed_load'} = res_bus.annual{:, 'uncurtailed_load'} - res_bus.annual{:, 'load'};
res_bus.annual{:, 'peak'} = max(load_hourly, [], 2);
res_bus.annual{:, 'lmp'} = sum((load_hourly.*hours').*lmp_hourly, 2) ./ res_bus.annual{:, 'load'};
idx = res_bus.annual{:, 'load'} == 0; % Re-compute zero loads lmp without load weighting
res_bus.annual{idx, 'lmp'} = lmp_hourly(idx, :) * probability;
res_bus.annual{:, 'lmp_prob'} = lmp_hourly(:, :) * probability;
%res_bus.annual{:, 'retail_prc'} = res_bus.annual{:, 'lmp'} * (1 + esc.dist_loss) + esc.dist_cost;
res_bus.annual{:, 'load_losses'} = res_bus.annual{:, 'load'} * esc.dist_loss;
res_bus.annual{:, 'retail_sales'} = res_bus.annual{:, 'load'} * (1 - esc.dist_loss);

%% Locational Information
if isfield(esc, 'bus_map')
    %loc_map = getMap(mpc, esc, 'bus', []); loc_map(:, 'bus') = [];
    %res_bus.annual = [res_bus.annual loc_map];
    res_bus.annual = join(res_bus.annual, esc.bus_map);
end

