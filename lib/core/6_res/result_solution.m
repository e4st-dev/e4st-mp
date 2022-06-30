function slv_res = result_solution(res)
% result_solution gathers information about solving and solution of LP

% E4ST
% Copyright (c) 2012-2022 by Power System Engineering Research Center (PSERC)
%
% This file is part of E4ST.
% Covered by the 3-clause BSD License (see LICENSE file for details).

slv_res.grb_out = table({num2str(res.success)}, 'VariableNames', {'Success'});
slv_res.grb_out{:, 'Status'} = {res.slv_status};
slv_res.grb_out{:, 'IterCount'} = res.opf_results.raw.output.itercount;
slv_res.grb_out{:, 'BarIterCount'} = res.opf_results.raw.output.baritercount;
%slv_res.grb_out{:, 'ObjVal'} = res.opf_results.raw.output.objval;

try
    slv_res.grb_out{:, 'ObjVal'} = res.opf_results.raw.output.objval;
catch
    slv_res.grb_out{:, 'ObjVal'} = -1;
end

slv_res.grb_out{:, {'check_toc_min', 'check_toc_max', 'check_caplim_min', ...
    'check_caplim_max', 'check_iflim_min', 'check_iflim_max'}} = 0;
if isfield(res, 'user_result') && isfield(res.user_result, 'constraints_res')
    if isfield(res.user_result.constraints_res, 'toc')
        toc = res.user_result.constraints_res.toc;
        if ~isempty(toc)
            slv_res.grb_out{:, 'check_toc_min'} = sum(toc{isfinite(toc{:, 'check_min'}), 'check_min'});
            slv_res.grb_out{:, 'check_toc_max'} = sum(toc{isfinite(toc{:, 'check_max'}), 'check_max'});
        end
    end
    if isfield(res.user_result.constraints_res, 'caplim')
        caplim = res.user_result.constraints_res.caplim;
        slv_res.grb_out{:, 'check_caplim_min'} = sum(caplim{isfinite(caplim{:, 'check_min'}), 'check_min'});
        slv_res.grb_out{:, 'check_caplim_max'} = sum(caplim{isfinite(caplim{:, 'check_max'}), 'check_max'});
    end
    if isfield(res.user_result.constraints_res, 'iflim')
        iflim = res.user_result.constraints_res.iflim;
        slv_res.grb_out{:, 'check_iflim_min'} = sum(iflim{isfinite(iflim{:, 'check_min'}), 'check_min'});
        slv_res.grb_out{:, 'check_iflim_max'} = sum(iflim{isfinite(iflim{:, 'check_max'}), 'check_max'});
    end
end
