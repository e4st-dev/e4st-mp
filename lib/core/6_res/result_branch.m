function res_branch = result_branch(res, mpc, esc)
%result_branch prepares branch-specific results  

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
% by Biao Mao (Rensselaer Polytechnic Institute) and Paul Picciano (Resources for the Future)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

%% Hourly Results
define_constants;
c = idx_dcline;

hours = esc.hrs_map{:, 'hours'};
probability = esc.hrs_map{:, 'probability'};
nh = length(hours);

% Allowing for cases with no DC lines.
if ~isfield(res.base, 'dcline')
    for i = 1:nh
        if i == 1
            res.base.dcline = zeros(0,23);
        else
            res.cont(i - 1).dcline = zeros(0,23);
        end
    end
end

flow_ac_hourly = zeros(size(res.base.branch, 1), nh);
flow_dc_hourly = zeros(size(res.base.dcline, 1), nh);
for i = 1:nh
    if i == 1
        branch = res.base.branch;
        dcline = res.base.dcline;
    else
        branch = res.cont(i - 1).branch;
        dcline = res.cont(i - 1).dcline;
    end
    flow_ac_hourly(:, i) = branch(:, PF); % PT = 16
    flow_dc_hourly(:, i) = dcline(:, c.PF); % PT = 16
    mu_ac_hourly(:, i) = branch(:, MU_SF) / probability(i);
    mu_min_dc_hourly(:, i) = dcline(:, c.MU_PMIN) / probability(i);
    mu_max_dc_hourly(:, i) = dcline(:, c.MU_PMAX) / probability(i);
end

% Order of index is exactly the same with mpc.bus
res_branch.flow_ac_hourly = flow_ac_hourly;
res_branch.flow_dc_hourly = flow_dc_hourly;
res_branch.mu_ac_hourly = mu_ac_hourly;
res_branch.mu_min_dc_hourly = mu_min_dc_hourly;
res_branch.mu_max_dc_hourly = mu_max_dc_hourly;

%% Annual Results
res_branch.annual_ac = array2table(mpc.branch(:, 1:2), 'VariableNames', {'fbus', 'tbus'});
res_branch.annual_ac{:, 'br_type'} = {'ac'};
res_branch.annual_ac{:, 'status'} = mpc.branch(:, BR_STATUS);
res_branch.annual_ac{:, 'RATE_A'} = mpc.branch(:, RATE_A);
res_branch.annual_ac{:, 'RATE_B'} = mpc.branch(:, RATE_B);
res_branch.annual_ac{:, 'RATE_C'} = mpc.branch(:, RATE_C);
% Add bus voltages
res_branch.annual_ac{:, 'flow_mwh'} = res_branch.flow_ac_hourly * hours;
res_branch.annual_ac{:, 'maxflow_mw'} = max(res_branch.flow_ac_hourly, [], 2);
res_branch.annual_ac{:, 'f_mu'} = sum((flow_ac_hourly.*hours').*mu_ac_hourly, 2) ./ res_branch.annual_ac{:, 'flow_mwh'};
idx = res_branch.annual_ac{:, 'flow_mwh'} == 0; % Re-compute zero loads lmp without load weighting
res_branch.annual_ac{idx, 'f_mu'} = mu_ac_hourly(idx, :) * probability; % Hour-weighting

% DC line annual results
if isfield(mpc, 'dcline')
    res_branch.annual_dc = array2table(mpc.dcline(:, 1:2), 'VariableNames', {'fbus', 'tbus'});
    res_branch.annual_dc{:, 'br_type'} = {'dc'};
    res_branch.annual_dc{:, 'flow_mwh'} = res_branch.flow_dc_hourly * hours;
    res_branch.annual_dc{:, 'maxflow_mw'} = max(res_branch.flow_dc_hourly, [], 2);
    res_branch.annual_dc{:, 'mu_min'} = sum((flow_dc_hourly.*hours').*mu_min_dc_hourly, 2) ./ res_branch.annual_dc{:, 'flow_mwh'};
    res_branch.annual_dc{:, 'mu_max'} = sum((flow_dc_hourly.*hours').*mu_max_dc_hourly, 2) ./ res_branch.annual_dc{:, 'flow_mwh'};
    idx = res_branch.annual_dc{:, 'flow_mwh'} == 0; % Re-compute zero loads lmp without load weighting
    res_branch.annual_dc{idx, 'mu_min'} = mu_min_dc_hourly(idx, :) * probability;
    res_branch.annual_dc{idx, 'mu_max'} = mu_max_dc_hourly(idx, :) * probability;
end

%% Merchandising Surplus
bus_res = res.user_result.bus_res;
% AC
flmp = join(res_branch.annual_ac(:, 'fbus'), [array2table(bus_res.annual{:, 'bus'}, 'VariableNames', {'fbus'}), array2table(bus_res.lmp_hourly)]);
tlmp = join(res_branch.annual_ac(:, 'tbus'), [array2table(bus_res.annual{:, 'bus'}, 'VariableNames', {'tbus'}), array2table(bus_res.lmp_hourly)]);
res_branch.ms_ac_hourly = (tlmp{:, 2:end} - flmp{:, 2:end}) .* res_branch.flow_ac_hourly;
res_branch.annual_ac{:, 'merch_surplus'} = res_branch.ms_ac_hourly * hours;

% DC
if isfield(mpc, 'dcline')
    flmp = join(res_branch.annual_dc(:, 'fbus'), [array2table(bus_res.annual{:, 'bus'}, 'VariableNames', {'fbus'}), array2table(bus_res.lmp_hourly)]);
    tlmp = join(res_branch.annual_dc(:, 'tbus'), [array2table(bus_res.annual{:, 'bus'}, 'VariableNames', {'tbus'}), array2table(bus_res.lmp_hourly)]);
    res_branch.ms_dc_hourly = (tlmp{:, 2:end} - flmp{:, 2:end}) .* res_branch.flow_dc_hourly;
    res_branch.annual_dc{:, 'merch_surplus'} = res_branch.ms_dc_hourly * hours;
end
