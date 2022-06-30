function offer = e4st_offer2mat(offer0, ng)
%E4ST_OFFER2MAT  Converts the various offer inputs into an offer matrix.
%   OFFER = E4ST_OFFER2MAT(OFFER0, NG)
%   OFFER = E4ST_OFFER2MAT(OFFER0)
%   Converts the offer data for c3sopf() or sopf2() from any of the
%   three input options (offer matrix, offer struct, mpc struct) into an
%   offer matrix. The first argument OFFER0 is an offer matrix, offer struct
%   or MPC struct containing offer data. The second argument NG is the
%   number of rows in the gen matrix and optional except when OFFER0 is an
%   offer struct. In this case, OFFER0 has the following fields for active
%   power quantities, each an ng x 1 vector ...
%       offer
%           .PositiveActiveReservePrice
%           .PositiveActiveReserveQuantity
%           .NegativeActiveReservePrice
%           .NegativeActiveReserveQuantity
%           .PositiveActiveDeltaPrice
%           .NegativeActiveDeltaPrice
%           .PositiveActiveReservePrice2        (optional quadratic term)
%           .NegativeActiveReservePrice2        (optional quadratic term)
%           .PositiveActiveDeltaPrice2          (optional quadratic term)
%           .NegativeActiveDeltaPrice2          (optional quadratic term)
%           .ActiveContractMin                  (optional)
%           .ActiveContractMax                  (optional)
%           .PminFactor                         (optional)
%   ... and optionally, the corresponding for reactive power ...
%           .PositiveReactiveReservePrice       (optional)
%           .PositiveReactiveReserveQuantity    (optional)
%           .NegativeReactiveReservePrice       (optional)
%           .NegativeReactiveReserveQuantity    (optional)
%           .PositiveReactiveDeltaPrice         (optional)
%           .NegativeReactiveDeltaPrice         (optional)
%           .PositiveReactiveReservePrice2      (optional quadratic term)
%           .NegativeReactiveReservePrice2      (optional quadratic term)
%           .PositiveReactiveDeltaPrice2        (optional quadratic term)
%           .NegativeReactiveDeltaPrice2        (optional quadratic term)
%           .ReactiveContractMin                (optional)
%           .ReactiveContractMax                (optional)
%   If OFFER0 is a matrix, the first NG rows contain the active power
%   quantities and the 2nd set of NG rows (optional) contain the reactive
%   power quantities. The columns correspond to the fields listed above
%   in the listed order. In this case, the matrix is simply returned as-is
%   with the exception that default values are added for any missing
%   optional columns.
%
%   If the OFFER0 argument is a MATPOWER case struct it should have
%   'reserve', 'energy_delta_cost' and 'contract' fields, which
%   take the following form, where offerp refers to the first NG rows of
%   the corresponding offer matrix and offerq to the optional 2nd set of
%   NG rows:
%       .reserve
%           .cost
%               .Rp_pos     [ offerp(:, 1) ]
%               .Rp_neg     [ offerp(:, 3) ]
%               .Rp_pos2    [ offerp(:, 7) ]    (optional quadratic term)
%               .Rp_neg2    [ offerp(:, 8) ]    (optional quadratic term)
%               .Rq_pos     [ offerq(:, 1) ]    (optional)
%               .Rq_neg     [ offerq(:, 3) ]    (optional)
%               .Rq_pos2    [ offerq(:, 7) ]    (optional quadratic term)
%               .Rq_neg2    [ offerq(:, 8) ]    (optional quadratic term)
%           .cap
%               .Rp_pos     [ offerp(:, 2) ]
%               .Rp_neg     [ offerp(:, 4) ]
%               .Rq_pos     [ offerq(:, 2) ]    (optional)
%               .Rq_neg     [ offerq(:, 4) ]    (optional)
%       .energy_delta_cost
%           .dP_pos         [ offerp(:, 5) ]
%           .dP_neg         [ offerp(:, 6) ]
%           .dP_pos2        [ offerp(:, 9) ]    (optional quadratic term)
%           .dP_neg2        [ offerp(:, 10)]    (optional quadratic term)
%           .dQ_pos         [ offerq(:, 5) ]    (optional)
%           .dQ_neg         [ offerq(:, 6) ]    (optional)
%           .dQ_pos2        [ offerq(:, 9) ]    (optional quadratic term)
%           .dQ_neg2        [ offerq(:, 10)]    (optional quadratic term)
%       .contract                               (optional)
%           .Pc_min         [ offerp(:, 11)]    (optional)
%           .Pc_max         [ offerp(:, 12)]    (optional)
%           .Qc_min         [ offerq(:, 11)]    (optional)
%           .Qc_max         [ offerq(:, 12)]    (optional)
%       .pmin_factor        [ offerp(:, 13)]    (optional)
%
%   The returned offer matrix has 12 columns unless OFFER0 is (1) a matrix
%   with 13 columns, (2) a case struct with a 'pmin_factor' field, or
%   (3) an offer struct with a 'PminFactor' field, in which case it has
%   13 columns.

%   E4ST
%   Copyright (c) 2000-2022 by Power System Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

if isstruct(offer0)
    if isfield(offer0, 'baseMVA')   %% offer0 is an MPC struct
        mpc = offer0;
        ng = size(mpc.gen, 1);
        nc = 12;
        if isfield(mpc, 'pmin_factor')
            nc = nc + 1;
        end
        offer = zeros(ng, nc);
        if (isfield(mpc, 'reserve') && ...  %% HAVE_Q ?
                ((isfield(mpc.reserve, 'cost') && ...
                    (isfield(mpc.reserve.cost, 'Rq_pos') || ...
                     isfield(mpc.reserve.cost, 'Rq_neg') || ...
                     isfield(mpc.reserve.cost, 'Rq_pos2') || ...
                     isfield(mpc.reserve.cost, 'Rq_neg2'))) || ...
                 (isfield(mpc.reserve, 'cap') && ...
                    (isfield(mpc.reserve.cap, 'Rq_pos') || ...
                     isfield(mpc.reserve.cap, 'Rq_neg'))))) || ...
            (isfield(mpc, 'energy_delta_cost') && ...
                (isfield(mpc.energy_delta_cost, 'dQ_pos') || ...
                 isfield(mpc.energy_delta_cost, 'dQ_neg') || ...
                 isfield(mpc.energy_delta_cost, 'dQ_pos2') || ...
                 isfield(mpc.energy_delta_cost, 'dQ_neg2'))) || ...
            (isfield(mpc, 'contract') && ...
                (isfield(mpc.contract, 'Qc_min') || ...
                 isfield(mpc.contract, 'Qc_max')))
          offerq = zeros(ng, 12);
        else
          offerq = [];
        end
        if isfield(mpc, 'reserve')
            if isfield(mpc.reserve, 'cost')
                if isfield(mpc.reserve.cost, 'Rp_pos')
                    offer(:, 1) = mpc.reserve.cost.Rp_pos;
                end
                if isfield(mpc.reserve.cost, 'Rp_neg')
                    offer(:, 3) = mpc.reserve.cost.Rp_neg;
                end
                if isfield(mpc.reserve.cost, 'Rp_pos2')
                    offer(:, 7) = mpc.reserve.cost.Rp_pos2;
                end
                if isfield(mpc.reserve.cost, 'Rp_neg2')
                    offer(:, 8) = mpc.reserve.cost.Rp_neg2;
                end
                if isfield(mpc.reserve.cost, 'Rq_pos')
                    offerq(:, 1) = mpc.reserve.cost.Rq_pos;
                end
                if isfield(mpc.reserve.cost, 'Rq_neg')
                    offerq(:, 3) = mpc.reserve.cost.Rq_neg;
                end
                if isfield(mpc.reserve.cost, 'Rq_pos2')
                    offerq(:, 7) = mpc.reserve.cost.Rq_pos2;
                end
                if isfield(mpc.reserve.cost, 'Rq_neg2')
                    offerq(:, 8) = mpc.reserve.cost.Rq_neg2;
                end
            end
            if isfield(mpc.reserve, 'cap')
                if isfield(mpc.reserve.cap, 'Rp_pos')
                    offer(:, 2) = mpc.reserve.cap.Rp_pos;
                end
                if isfield(mpc.reserve.cap, 'Rp_neg')
                    offer(:, 4) = mpc.reserve.cap.Rp_neg;
                end
                if isfield(mpc.reserve.cap, 'Rq_pos')
                    offerq(:, 2) = mpc.reserve.cap.Rq_pos;
                end
                if isfield(mpc.reserve.cap, 'Rq_neg')
                    offerq(:, 4) = mpc.reserve.cap.Rq_neg;
                end
            end
        end
        if isfield(mpc, 'energy_delta_cost')
            if isfield(mpc.energy_delta_cost, 'dP_pos')
                offer(:, 5) = mpc.energy_delta_cost.dP_pos;
            end
            if isfield(mpc.energy_delta_cost, 'dP_neg')
                offer(:, 6) = mpc.energy_delta_cost.dP_neg;
            end
            if isfield(mpc.energy_delta_cost, 'dP_pos2')
                offer(:, 9) = mpc.energy_delta_cost.dP_pos2;
            end
            if isfield(mpc.energy_delta_cost, 'dP_neg2')
                offer(:, 10) = mpc.energy_delta_cost.dP_neg2;
            end
            if isfield(mpc.energy_delta_cost, 'dQ_pos')
                offerq(:, 5) = mpc.energy_delta_cost.dQ_pos;
            end
            if isfield(mpc.energy_delta_cost, 'dQ_neg')
                offerq(:, 6) = mpc.energy_delta_cost.dQ_neg;
            end
            if isfield(mpc.energy_delta_cost, 'dQ_pos2')
                offerq(:, 9) = mpc.energy_delta_cost.dQ_pos2;
            end
            if isfield(mpc.energy_delta_cost, 'dQ_neg2')
                offerq(:, 10) = mpc.energy_delta_cost.dQ_neg2;
            end
        end
        if isfield(mpc, 'contract')
            if isfield(mpc.contract, 'Pc_min')
                offer(:, 11) = mpc.contract.Pc_min;
            else
                offer(:, 11) = -Inf;
            end
            if isfield(mpc.contract, 'Pc_max')
                offer(:, 12) = mpc.contract.Pc_max;
            else
                offer(:, 12) = Inf;
            end
            if isfield(mpc.contract, 'Qc_min')
                offerq(:, 11) = mpc.contract.Qc_min;
            else
                offerq(:, 11) = -Inf;
            end
            if isfield(mpc.contract, 'Qc_max')
                offerq(:, 12) = mpc.contract.Qc_max;
            else
                offerq(:, 12) =  Inf;
            end
        end
        if isfield(mpc, 'pmin_factor')
            offer(:, 13) = mpc.pmin_factor;
            if ~isempty(offerq)
                offerq(:, 13) = 0;
            end
        end
        if ~isempty(offerq)
            offer = [ offer; offerq ];
        end
    else                            %% offer0 is an offer struct
        if nargin < 2
            error('e4st_offer2mat: requires 2nd arg (ng) when called with offer struct');
        end
        nc = 12;
        if isfield(offer0, 'PminFactor')
            nc = nc + 1;
        end
        offer = zeros(ng, nc);
        if isfield(offer0, 'PositiveActiveReservePrice')
            offer(:, 1) = offer0.PositiveActiveReservePrice;
        end
        if isfield(offer0, 'PositiveActiveReserveQuantity')
            offer(:, 2) = offer0.PositiveActiveReserveQuantity;
        end
        if isfield(offer0, 'NegativeActiveReservePrice')
            offer(:, 3) = offer0.NegativeActiveReservePrice;
        end
        if isfield(offer0, 'NegativeActiveReserveQuantity')
            offer(:, 4) = offer0.NegativeActiveReserveQuantity;
        end
        if isfield(offer0, 'PositiveActiveDeltaPrice')
            offer(:, 5) = offer0.PositiveActiveDeltaPrice;
        end
        if isfield(offer0, 'NegativeActiveDeltaPrice')
            offer(:, 6) = offer0.NegativeActiveDeltaPrice;
        end
        if isfield(offer0, 'PositiveActiveReservePrice2')
            offer(:, 7) = offer0.PositiveActiveReservePrice2;
        end
        if isfield(offer0, 'NegativeActiveReservePrice2')
            offer(:, 8) = offer0.NegativeActiveReservePrice2;
        end
        if isfield(offer0, 'PositiveActiveDeltaPrice2')
            offer(:, 9) = offer0.PositiveActiveDeltaPrice2;
        end
        if isfield(offer0, 'NegativeActiveDeltaPrice2')
            offer(:, 10) = offer0.NegativeActiveDeltaPrice2;
        end
        if isfield(offer0, 'ActiveContractMin')
            offer(:, 11) = offer0.ActiveContractMin;
        else
            offer(:, 11) = -Inf;
        end
        if isfield(offer0, 'ActiveContractMax')
            offer(:, 12) = offer0.ActiveContractMax;
        else
            offer(:, 12) =  Inf;
        end
        if isfield(offer0, 'PminFactor')
            offer(:, 13) = offer0.PminFactor;
        end
        if isfield(offer0, 'PositiveReactiveReservePrice') || ...
                isfield(offer0, 'PositiveReactiveReserveQuantity') || ...
                isfield(offer0, 'NegativeReactiveReservePrice') || ...
                isfield(offer0, 'NegativeReactiveReserveQuantity') || ...
                isfield(offer0, 'PositiveReactiveDeltaPrice') || ...
                isfield(offer0, 'NegativeReactiveDeltaPrice') || ...
                isfield(offer0, 'PositiveReactiveReservePrice2') || ...
                isfield(offer0, 'NegativeReactiveReservePrice2') || ...
                isfield(offer0, 'PositiveReactiveDeltaPrice2') || ...
                isfield(offer0, 'NegativeReactiveDeltaPrice2') || ...
                isfield(offer0, 'ReactiveContractMin') || ...
                isfield(offer0, 'ReactiveContractMax')
            offerq = zeros(ng, 12);
            if isfield(offer0, 'PositiveReactiveReservePrice')
                offerq(:, 1) = offer0.PositiveReactiveReservePrice;
            end
            if isfield(offer0, 'PositiveReactiveReserveQuantity')
                offerq(:, 2) = offer0.PositiveReactiveReserveQuantity;
            end
            if isfield(offer0, 'NegativeReactiveReservePrice')
                offerq(:, 3) = offer0.NegativeReactiveReservePrice;
            end
            if isfield(offer0, 'NegativeReactiveReserveQuantity')
                offerq(:, 4) = offer0.NegativeReactiveReserveQuantity;
            end
            if isfield(offer0, 'PositiveReactiveDeltaPrice')
                offerq(:, 5) = offer0.PositiveReactiveDeltaPrice;
            end
            if isfield(offer0, 'NegativeReactiveDeltaPrice')
                offerq(:, 6) = offer0.NegativeReactiveDeltaPrice;
            end
            if isfield(offer0, 'PositiveReactiveReservePrice2')
                offerq(:, 7) = offer0.PositiveReactiveReservePrice2;
            end
            if isfield(offer0, 'NegativeReactiveReservePrice2')
                offerq(:, 8) = offer0.NegativeReactiveReservePrice2;
            end
            if isfield(offer0, 'PositiveReactiveDeltaPrice2')
                offerq(:, 9) = offer0.PositiveReactiveDeltaPrice2;
            end
            if isfield(offer0, 'NegativeReactiveDeltaPrice2')
                offerq(:, 10) = offer0.NegativeReactiveDeltaPrice2;
            end
            if isfield(offer0, 'ReactiveContractMin')
                offerq(:, 11) = offer0.ReactiveContractMin;
            else
                offerq(:, 11) = -Inf;
            end
            if isfield(offer0, 'ReactiveContractMax')
                offerq(:, 12) = offer0.ReactiveContractMax;
            else
                offerq(:, 12) =  Inf;
            end
            if isfield(offer0, 'PminFactor')
                offerq(:, 13) = 0;
            end
            offer = [ offer; offerq ];
        end
    end
else                                %% offer0 is a matrix
    offer = offer0;
    [m, n] = size(offer);
    if n < 12       %% pad on right if necessary
        offer = [offer zeros(m, 12-n)];
        offer(:, 12) = Inf;
        if n < 11
            offer(:, 11) = -Inf;
        end
    end
end
