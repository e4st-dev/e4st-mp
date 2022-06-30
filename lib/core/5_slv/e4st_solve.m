function [results, f, success, info, et, g, jac, xr, pimul] = e4st_solve(varargin)
%E4ST_SOLVE  Core solver for E4ST.
%
%   Examples:
%       Output argument options:
%
%       [results, f, success] = e4st_solve(...)
%       [results, f, success, info] = e4st_solve(...)
%       [results, f, success, info, et] = e4st_solve(...)
%       [results, f, success, info, et, g, jac] = e4st_solve(...)
%       [results, f, success, info, et, g, jac, xr, pimul] = e4st_solve(...)
%
%       Input argument options:
%
%       e4st_solve(mpc)
%       e4st_solve(mpc, mpopt)
%       e4st_solve(mpc, offer, contab)
%       e4st_solve(mpc, offer, contab, mpopt)
%
%   mpc is a MATPOWER case file or case struct with the fields baseMVA, bus,
%   gen, branch, and (optionally) areas. It may also include a 'contingencies'
%   field (in place of contab argument). The offer argument can be a struct
%   or a matrix. If it is a struct, it has the following fields for active
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
%   If offer is a matrix, the first ng rows contain the active power
%   quantities and the 2nd set of ng rows (optional) contain the reactive
%   power quantities. The columns correspond to the fields listed above
%   in the listed order.
%
%   Alternatively, the offer argument can be omitted and the fields
%   'reserve', 'energy_delta_cost' and 'contract' included in mpc.
%   In this case, the 'reserve', 'energy_delta_cost' and 'contract' fields
%   take the following form, where offerp refers to the first ng rows of
%   the corresponding offer matrix and offerq to the optional 2nd set of
%   ng rows:
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
%   An optional 'availability_factor' field can be used to specify an
%   availability factor for each generator in each scenario. It can be
%   either an ng x 1 vector or ng x (nc+1) matrix, where each element is
%   between 0 and 1.
%
%   contab is the contingency table, type 'help apply_changes' and
%   'help idx_ct' for details about the format.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.

%   E4ST
%   Copyright (c) 2000-2022 by Power System Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of E4ST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

% SECTION 1: ARGUMENT PARSING, REORDERING, UNPACKING

[baseMVA, bus, gen, branch, gencost, dcline, iflims, softlims, offer, ...
    contab, mpopt, HAVE_Q, caplim_map, caplim_max, caplim_min, avail_fac, ...
    toc_map, toc_max, toc_min, toc_coeff, toc_type, ...
    days, short_term_storage] = e4st_args(varargin{:});
N = [];

if nargout > 5
    mpopt = mpoption(mpopt, 'opf.return_raw_der', 1);
end

% options
OUT_ALL = mpopt.out.all;
dc      = strcmp(upper(mpopt.model), 'DC');
if isfield(mpopt, 'sopf') && isfield(mpopt.sopf, 'force_Pc_eq_P0')
    FORCE_PC_EQ_P0 = mpopt.sopf.force_Pc_eq_P0;
else
    FORCE_PC_EQ_P0 = 0;     %% off by default
end

if mpopt.verbose > 0
  v = e4st_ver('all');
  fprintf('\nE4ST Version %s, %s\n', v.Version, v.Date);
  fprintf('Engineering, Economic, and Environmental Electricity Simulation Tool\n');
end

if dc                   %% force HAVE_Q to false for DC runs
    HAVE_Q = 0;
end

% Load column indices for case tables.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;
c = idx_dcline;

% If tables do not have multiplier/extra columns, append zero cols.
% Update whenever the data format changes!
if size(bus,2) < MU_VMIN
  bus = [bus zeros(size(bus,1),MU_VMIN-size(bus,2)) ];
end
if size(gen,2) < MU_QMIN
  gen = [ gen zeros(size(gen,1),MU_QMIN-size(gen,2)) ];
end
if size(branch,2) < MU_ANGMAX
  branch = [ branch zeros(size(branch,1),MU_ANGMAX-size(branch,2)) ];
end
if ~isempty(dcline) && size(dcline,2) < c.MU_QMAXT
  dcline = [ dcline zeros(size(dcline,1),c.MU_QMAXT-size(dcline,2)) ];
end

% get number of contingencies & list of labels
[xx, ii] = sort(contab(:, CT_LABEL)); %sort in ascending contingency label
contab = contab(ii, :);
clist0 = unique(contab(:, CT_LABEL));
if isempty(days)
    HAVE_DAYS = 0;
    clist = clist0;
else
    HAVE_DAYS = 1;
    clist = vertcat(days{:, 1});
    if clist(1) == 0
        clist(1) = [];  %% remove "base case" from list
    else
        error('e4st_solve: first hour in first day in ''days'' field must be 0 (base case).')
    end
    if any(~ismember(clist, clist0))
        error('e4st_solve: ''days'' field contains invalid contingency labels');
    end
end

% Filter-out uncomitted equipment from data; the row indices in
% contab must be modified accordingly.  This might lead to
% deleting contingencies... yuk... must track this too if this
% is to be helpful to a higher-level decommitment routine...

% More on higher-level unit de-commitment: an outer de-commit
% decision for a generator or a line might make a corresponding
% contingency in the original list irrelevant; the problem becomes
% entirely different, with different base case probability and so on.
% To keep things consistent, on output, the post-contingency
% flows for deleted contingencies should still have a numbered spot
% in the data but be empty, which must be correctly interpreted
% by the calling unit de-commitment code.

base_on_gen  = find(gen(:, GEN_STATUS) >  0);
base_off_gen = find(gen(:, GEN_STATUS) <= 0);
base_on_branch  = find(branch(:, BR_STATUS) >  0);
base_off_branch = find(branch(:, BR_STATUS) <= 0);
gen_original = gen;                 % save these three
branch_original = branch;
clist_original = clist;
gen =  gen(base_on_gen, :);         % now these contain only equipment
branch = branch(base_on_branch, :); % committed on input
ng = size(gen_original, 1);         % original size(gen), at least temporarily
if size(gencost,1) == ng            % while we do something with it
  gencost = gencost(base_on_gen, :);
else
  gencost = gencost( [base_on_gen; base_on_gen+ng], :);
end
if ~isempty(short_term_storage)
  nsts0 = size(short_term_storage, 1);  %% original number of short-term storage units
  ge2i = zeros(ng, 1);      %% create external-to-internal gen index map
  ge2i(base_on_gen) = 1:length(base_on_gen);
  base_sts_on = find(ismember(short_term_storage(:, 1), base_on_gen));
  base_sts_off = find(ismember(short_term_storage(:, 1), base_off_gen));
  short_term_storage_original = short_term_storage;
  short_term_storage = short_term_storage(base_sts_on, :);
  short_term_storage(:, 1) = ge2i(short_term_storage(:, 1));
  istsp = short_term_storage(:, 1); % (internal) gen indices of short-term storage units
  xi = short_term_storage(:, 2);    % efficiencies of short-term storage units
  tsts = short_term_storage(:, 3);  % time to full discharge of short-term storage units
  nsts = length(istsp);             % number of short term storage units
  % error if a storage unit does not have negative PMIN and positive PMAX
  if any(gen(istsp, PMIN) >= 0 | gen(istsp, PMAX) <= 0)
      error('e4st_solve: short term storage units must have PMIN < 0 and PMAX > 0');
  end
  % error if a contingency turns off a short term storage unit
  k = find(contab(:, CT_TABLE) == CT_TGEN & ...
          ismember(contab(:, CT_ROW), [0; short_term_storage(:, 1)]) & ...
          contab(:, CT_COL) == GEN_STATUS & ...
          contab(:, CT_CHGTYPE) == CT_REP & ...
          contab(:, CT_NEWVAL) <= 0 );
  if ~isempty(k)
      error('e4st_solve: switching off short term storage units via CONTAB is not supported');
  end
else
  nsts = 0;
end
if ~isempty(dcline)
  base_on_dcline  = find(dcline(:, c.BR_STATUS) >  0);
  base_off_dcline = find(dcline(:, c.BR_STATUS) <= 0);
  dcline_original = dcline;
  dcline = dcline(base_on_dcline, :);
end
if ~isempty(iflims)
  ifmap = iflims.map;
  ifmap_original = ifmap;
  % update branch indices in ifmap, remove lines that are out-of-service
  e2i = zeros(size(branch_original, 1), 1);
  e2i(base_on_branch) = (1:size(branch, 1))';
  d = sign(ifmap(:, 2));
  br = abs(ifmap(:, 2));
  ifmap(:, 2) = d .* e2i(br);
  ifmap(ifmap(:, 2) == 0, :) = [];  %% delete branches that are out
end
if ~isempty(softlims)
    error('e4st_solve: softlims implementation has no tests and has not been updated for MATPOWER 7');

    softlimsmap = softlims.idx;
    softlimsmap_original = softlimsmap;
    % update branch indices in softlimsmap, remove lines that are out-of-service
    e2i = zeros(size(branch_original, 1), 1);
    e2i(base_on_branch) = (1:size(branch, 1))';
    d = sign(softlimsmap(:, 1));
    br = abs(softlimsmap(:, 1));
    softlimsmap(:, 1) = d .* e2i(br);
    softlimsmap(softlimsmap(:, 1) == 0, :) = []; %% delete branches that are out
end
if ~isempty(caplim_map)
  caplim_map = caplim_map(:, base_on_gen);
end
if ~isempty(avail_fac)
  avail_fac = avail_fac(base_on_gen, :);        % reduce to committed gens
  af_nc = size(avail_fac, 2);
  if af_nc > 1 && af_nc ~= length(clist) + 1    % check for proper # of cols
    error('e4st_solve.m: # of columns in mpc.availability_factor (%d) must equal 1 or # contingencies + 1 (%d)', ...
        af_nc, length(clist) + 1);
  end
end
if ~isempty(toc_map)
  toc_map = toc_map(:, base_on_gen, :);
  toc_coeff = toc_coeff(base_on_gen, :);
end

% Do the same to offer table/struct
% offer_original = offer;
if HAVE_Q
  offer = offer([base_on_gen; base_on_gen+ng], :);
else
  offer = offer(base_on_gen, :);
end

% Look for contingency table rows that act on branches or generators that
% are turned off on input (and therefore have now been removed). If so, delete
% the corresponding contingency table row. At the end, see if a contingency
% label in clist points to (now) nonexistent labels in the contingency
% table and if so entirely delete the contingency label from clist.
% NOTE: Does not catch changes that act on all rows using a row value of 0,
% or e.g. all gens in an area... yet ***************

rowcomlist = ones(size(contab,1), 1);
for l = 1:size(contab, 1)
  if contab(l, CT_ROW) ~= 0 && ...
        ((contab(l, CT_TABLE) == CT_TGEN && ...
            gen_original(contab(l, CT_ROW), GEN_STATUS) <= 0) || ...
         (contab(l, CT_TABLE) == CT_TBRCH && ...
            branch_original(contab(l, CT_ROW), BR_STATUS) <= 0))
    rowcomlist(l) = 0;
  end
end
contab = contab(rowcomlist ~= 0, :);
clabelcomlist = ones(size(clist));
for l = 1:length(clist)
  if ~any(clist(l) == contab(:, CT_LABEL))
    clabelcomlist(l) = 0;
  end
end
if any(clabelcomlist == 0)
    if HAVE_DAYS
        error('e4st_solve: does not support ''days'' with hours that modify off-line (eliminated) units.')
    else
        clist = clist(clabelcomlist ~= 0);
        % clabeldecomlist = clist(clabelcomlist == 0); % contingencies with these labels were thrown out
    end
end
nc = length(clist);

% resize avail_fac to make it ng x (nc+1)
if ~isempty(avail_fac)
  if size(avail_fac, 2) > 1
    avail_fac = avail_fac(:, [1; clabelcomlist] ~= 0);  % reduce for removed contingencies
  else
    avail_fac = avail_fac * ones(1, nc+1);  % expand vector to matrix
  end
else
  avail_fac = ones(ng, nc+1);               % all ones if not given
end

% Now catch renumbering of CT_ROW in contab due to deletion of rows
% for equipment that came with STATUS=off on input. Remember that at this
% point contab points only to rows that denote equipment that is to be kept
% in the analysis.

newrowgen = cumsum(gen_original(:,GEN_STATUS) > 0);   % newrowgen(4) contains
newrowbrch = cumsum(branch_original(:,BR_STATUS) > 0);% new row index for gen 4,
if size(gencost,1) == ng                              % assuming gen4 is kept
  newrowgencost = newrowgen;
else
  newrowgencost = cumsum([  gen_original(:,GEN_STATUS);
                            gen_original(:,GEN_STATUS) ]  > 0);
end
for l = 1:size(contab,1)
  if contab(l, CT_TABLE) == CT_TGEN && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowgen(contab(l, CT_ROW));
  elseif contab(l, CT_TABLE) == CT_TGENCOST && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowgencost(contab(l, CT_ROW));
  elseif contab(l, CT_TABLE) == CT_TBRCH && contab(l, CT_ROW) > 0
    contab(l, CT_ROW) = newrowbrch(contab(l, CT_ROW));
  end
end

% compute probability of nominal (base) flow as the complement of the
% sum of the probabilities of the contingencies
prob0 = 1;
for i=1:nc
  tmp = contab(clist(i) == contab(:, CT_LABEL), CT_PROB);
  prob0 = prob0 - tmp(1);
end

% Renumber buses consecutively (rows themselves in bus() are not permuted)
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
if ~isempty(dcline)
  e2i = sparse(max(i2e), 1);
  e2i(i2e) = (1:size(bus, 1))';
  dcline(:, c.F_BUS) = e2i( dcline(:, c.F_BUS) );
  dcline(:, c.T_BUS) = e2i( dcline(:, c.T_BUS) );
end

ng = size(gen,1);

% split short term storage into positive/negative generators
if nsts
    if HAVE_Q
        error('e4st_solve: support for short term storage with reactive power offers not implemented.');
    end

    % duplicate gens
    sts_gen = gen(istsp, :);    % negative part
    % PMIN of positive, PMAX of negative to zero
    gen(istsp, PMIN) = 0;
    sts_gen(:, PMAX) = eps;     % to avoid being handled as a dispatchable load
    istsm = (ng+1:ng+nsts)';    % indices of negative part
    gen = [gen; sts_gen];
    
    % duplicate gencosts
    sts_gencost = gencost(istsp, :);
    % zero out costs for negative part
    sts_gencost(:, MODEL) = POLYNOMIAL;
    sts_gencost(:, NCOST) = 1;
    sts_gencost(:, COST:end) = 0;
    gencost = [gencost; sts_gencost];

    % duplicate offers
    offer(istsp, [11 12]) = 0;      % force ActiveContractMin/Max to 0 for both parts
    sts_offer = offer(istsp, :);
    % zero out reserve price of negative part
    sts_offer(:, [1 3]) = 0;
    % set NegativeActiveReserveQuantity of negative part by PositiveActiveReserveQuantity
    sts_offer(:, 4) = offer(istsp, 2);
    % no need to worry about PminFactor for negative part, it is explicitly
    % excluded from PminFactor constraints
    offer = [offer; sts_offer];

    % expand other inputs
    if ~isempty(caplim_map) % negative part excluded from cap lims
      caplim_map = [caplim_map zeros(size(caplim_map, 1), nsts)];
    end
    if ~isempty(avail_fac)  % same availability factors for negative part
      avail_fac = [avail_fac; avail_fac(istsp, :)];
    end
    if ~isempty(toc_map)
      toc_map(:, istsm, :) = 0;     % negative part excluded from total output constraints
      toc_coeff = [toc_coeff; zeros(nsts, size(toc_coeff, 2))];
    end
    
    %% update ng
    ng = size(gen, 1);
end

[gencost_p, gencost_q] = pqcost(gencost, ng);

% Offer processing - from here on we have unified source of offer data
PositiveActiveReservePrice      = offer(1:ng, 1);
PositiveActiveReserveQuantity   = offer(1:ng, 2);
NegativeActiveReservePrice      = offer(1:ng, 3);
NegativeActiveReserveQuantity   = offer(1:ng, 4);
PositiveActiveDeltaPrice        = offer(1:ng, 5);
NegativeActiveDeltaPrice        = offer(1:ng, 6);
PositiveActiveReservePrice2     = offer(1:ng, 7);
NegativeActiveReservePrice2     = offer(1:ng, 8);
PositiveActiveDeltaPrice2       = offer(1:ng, 9);
NegativeActiveDeltaPrice2       = offer(1:ng, 10);
ActiveContractMin               = offer(1:ng, 11);
ActiveContractMax               = offer(1:ng, 12);
if size(offer, 2) >= 13
    PminFactor                  = offer(1:ng, 13);
else
    PminFactor                  = zeros(ng, 1);
end
if HAVE_Q
  PositiveReactiveReservePrice      = offer(ng+1:2*ng, 1);
  PositiveReactiveReserveQuantity   = offer(ng+1:2*ng, 2);
  NegativeReactiveReservePrice      = offer(ng+1:2*ng, 3);
  NegativeReactiveReserveQuantity   = offer(ng+1:2*ng, 4);
  PositiveReactiveDeltaPrice        = offer(ng+1:2*ng, 5);
  NegativeReactiveDeltaPrice        = offer(ng+1:2*ng, 6);
  PositiveReactiveReservePrice2     = offer(ng+1:2*ng, 7);
  NegativeReactiveReservePrice2     = offer(ng+1:2*ng, 8);
  PositiveReactiveDeltaPrice2       = offer(ng+1:2*ng, 9);
  NegativeReactiveDeltaPrice2       = offer(ng+1:2*ng, 10);
  ReactiveContractMin               = offer(ng+1:2*ng, 11);
  ReactiveContractMax               = offer(ng+1:2*ng, 12);
end
if any(any(offer(:, 7:10))) % we have quadratic inc/dec/reserve cost terms
  HAVE_QUADRATIC = 1;
else
  HAVE_QUADRATIC = 0;
end
if FORCE_PC_EQ_P0 && any(any(isfinite(offer(:, 11:12))))
  fprintf('\ne4st_solve: WARNING: Using limits on Pc (or Qc) is not recommended\n');
  fprintf('                 when the ''sopf.force_Pc_eq_P0'' option is enabled.\n\n');
end

% SECTION 2: CONTINGENCY PROCESSING AND SYSTEM AUGMENTATION

% The big-system building must be a two-pass process since
% start/end indices of both variables and constraints cannot be known
% until all contingencies are processed and we see how many
% injections and branches are active in each flow.  After the first
% pass, indices can be computed and then additional linear constraints
% can be built in the second pass.

% FIRST PASS. Augmented network information is built, as well as
% basic augmented generator and cost tables.
ng(1) = size(gen, 1);
nb(1) = size(bus, 1);   % not supposed to change, but..
nl(1) = size(branch, 1);
ndc(1) = size(dcline, 1);
if ~isempty(iflims)
    nif(1) = size(iflims.lims, 1);
    nifm(1) = size(ifmap, 1);
    ifidmax(1) = max(ifmap(:, 1));
else
    nif(1) = 0;
    nifm(1) = 0;
    ifidmax(1) = 0;
end
if ~isempty(softlims)
    nsoftlims(1) = size(softlims.cost, 1);
    nsoftlimsm(1) = size(softlimsmap, 1);
    softlimsidmax(1) = length(softlimsmap(:, 1));
else
    nsoftlims(1) = 0;
    nsoftlimsm(1) = 0;
    softlimsidmax(1) = 0;
end
Augmented_bus = bus;
Augmented_branch = branch;
Augmented_gen = gen;
Augmented_gencost_p = modcost(gencost_p, prob0);
if ~isempty(gencost_q)
  Augmented_gencost_q = modcost(gencost_q, prob0);
else
  Augmented_gencost_q = [];
end
if ~isempty(dcline)
  Augmented_dcline = dcline;
end
if ~isempty(iflims)
    Augmented_ifmap = ifmap;
    Augmented_iflims = iflims.lims(:, :, 1);
end
if ~isempty(softlims)
    Augmented_softlimsmap = softlimsmap;
    Augmented_softlimscost = softlims.cost(:, :, 1);
end
gen_stat = [];
branch_stat = [];
% dcline_stat = [];

if ~HAVE_DAYS
  days = num2cell([0; clist]);
end
ndays = size(days, 1);
k = 0;
for d = 1:ndays
 hrs = days{d, 1};
 if d == 1
  h0 = 2;
 else
  h0 = 1;
 end
 for h = h0:length(hrs)
  k = k + 1;
  clabel = clist(k);
  kthcontab = contab(clabel == contab(:, CT_LABEL), :); % get rows addressing kth contingency
  ng(k+1) = ng(1); % initialize to 'full size' and we'll decrement from there
  nb(k+1) = nb(1);
  nl(k+1) = nl(1);
  ndc(k+1) = ndc(1);
  nif(k+1) = nif(1);
  nifm(k+1) = nifm(1);
  ifidmax(k+1) = ifidmax(k) + ifidmax(1);
  softlimsidmax(k+1) = softlimsidmax(k) + softlimsidmax(1);

  % apply the modifications
  mpck = struct('bus', bus, 'gen', gen, 'branch', branch, 'gencost', gencost);
  mpck = apply_changes(clabel, mpck, contab);

  % working copies of tables for k-th contingency
  [buswork, genwork, branchwork] = deal(mpck.bus, mpck.gen, mpck.branch);
  [gencost_p_work, gencost_q_work] = pqcost(mpck.gencost, ng(1));

  % now catch any dimension changes that result from them and construct
  % augmented network data for this contingency
  nb(k+1) = size(bus,1);
  busbastmp = sum(nb(1:k));  % total # buses in base and previous flows
  buswork(:, BUS_I) = bus(:, BUS_I) + busbastmp;
  % augmented problem has multiple REF buses
  % opf() has been updated to fix all such REF bus voltage angles and
  % emit a warning in case that was not intentional
  gen_stat = [ gen_stat genwork(:, GEN_STATUS) ];
  ii = find(genwork(:, GEN_STATUS) <= 0);
  genwork(ii, : ) = [];
  gencost_p_work(ii, :) = [];
  if ~isempty(gencost_q_work)
    gencost_q_work(ii, : ) = [];
  end
  ng(k+1) = size(genwork, 1);
  genwork(:, GEN_BUS) = genwork(:, GEN_BUS) + busbastmp; % bus #s kth island
  branch_stat = [ branch_stat branchwork(:,BR_STATUS) ];
%  branchwork(:, RATE_A) = branchwork(:, RATE_B);
  branchwork(branchwork(:, BR_STATUS) <= 0, : ) = [];
  nl(k+1) = size(branchwork, 1);
  brbastmp = sum(nl(1:k));  % total # branches in base and previous flows
  branchwork(:, F_BUS) = branchwork(:, F_BUS) + busbastmp;
  branchwork(:, T_BUS) = branchwork(:, T_BUS) + busbastmp;
  if ~isempty(dcline)
    dclinework = dcline;    %% assume no changes per contingency (for now)
%     dcline_stat = [ dcline_stat dclinework(:,BR_STATUS) ];
%     dclinework(dclinework(:, BR_STATUS) <= 0, : ) = [];
%     ndc(k+1) = size(dclinework, 1);
    dclinework(:, c.F_BUS) = dclinework(:, c.F_BUS) + busbastmp;
    dclinework(:, c.T_BUS) = dclinework(:, c.T_BUS) + busbastmp;
    Augmented_dcline = [Augmented_dcline; dclinework];
  end
  if ~isempty(iflims)
    ifmapwork = ifmap;          %% assume no changes per contingency
    iflimswork = iflims.lims(:, :, k+1); %% assume there ARE changes per contingency
    ifbastmp = ifidmax(k);      %% max if id num in previous flows
    ifmapwork(:, 1)  = ifmapwork(:, 1) + ifbastmp;
    iflimswork(:, 1) = iflimswork(:, 1) + ifbastmp;
    e2i = zeros(size(mpck.branch, 1), 1);
    e2i(mpck.branch(:, BR_STATUS) >  0) = (1:size(branchwork, 1))';
    d = sign(ifmapwork(:, 2));
    br = abs(ifmapwork(:, 2));
    ifmapwork(:, 2) = d .* (e2i(br) + brbastmp);
    ifmapwork(ifmapwork(:, 2) == brbastmp, :) = [];  %% delete branches that are out
    nifm(k+1) = size(ifmapwork, 1);
    Augmented_ifmap = [Augmented_ifmap; ifmapwork];
    Augmented_iflims = [Augmented_iflims; iflimswork];
  end
  if ~isempty(softlims)
      softlimsmapwork = softlimsmap; %% assume no changes per contingency
      softlimscostwork = softlims.cost(:, :, k+1); %% assume no changes per contingency
      %softlimsbastmp = softlimsidmax(k);      %% max if id num in previous flows
      %softlimsmapwork(:, 1)  = softlimsmapwork(:, 1) + softlimsbastmp;
      %softlimscostwork(:, 1) = softlimscostwork(:, 1) + softlimsbastmp;
      e2i = zeros(size(mpck.branch, 1), 1);
      e2i(mpck.branch(:, BR_STATUS) > 0) = (1:size(branchwork, 1))';
      d = sign(softlimsmapwork(:, 1));
      br = abs(softlimsmapwork(:, 1));
      softlimsmapwork(:, 1) = d .* (e2i(br) + brbastmp);
      softlimsmapwork(softlimsmapwork(:, 1) == brbastmp, :) = []; %% delete branches that are out
      nsoftlimsm(k+1) = size(softlimsmapwork, 1);
      Augmented_softlimsmap = [Augmented_softlimsmap; softlimsmapwork];
      Augmented_softlimscost = [Augmented_softlimscost; softlimscostwork];
  end
  Augmented_bus = [Augmented_bus; buswork];
  Augmented_gen = [Augmented_gen; genwork];
  Augmented_branch = [Augmented_branch; branchwork];
  Augmented_gencost_p = [Augmented_gencost_p; ...
                         modcost(gencost_p_work, kthcontab(1, CT_PROB)) ];
  Augmented_gencost_q = [Augmented_gencost_q; ...
                         modcost(gencost_q_work, kthcontab(1, CT_PROB)) ];
 end    % for h=hours
end     % for d=days


Augmented_gencost = [ Augmented_gencost_p; Augmented_gencost_q];

% SECTION 3: INDEXING SCHEME FOR VARIABLES AND NONLINEAR CONSTRAINTS

% build variable indexing scheme
nb_total = sum(nb(:));
ng_total = sum(ng(:));
nl_total = sum(nl(:));
% ndc_total = sum(ndc(:));
thbas(1) = 1;                      thend(1) = thbas(1) + nb(1) - 1;
vbas(1)  = thbas(1) + nb_total;    vend(1)  = vbas(1)  + nb(1) - 1;
pgbas(1) = 2*nb_total + 1;         pgend(1) = pgbas(1) + ng(1) - 1;
qgbas(1) = 2*nb_total+ng_total+1;  qgend(1) = qgbas(1) + ng(1) - 1;
if nc > 0
  for k = 2:nc+1
    thbas(k) = thend(k-1) + 1;     thend(k) = thbas(k) + nb(k) - 1;
    vbas(k)  = vend(k-1)  + 1;     vend(k)  = vbas(k)  + nb(k) - 1;
    pgbas(k) = pgend(k-1) + 1;     pgend(k) = pgbas(k) + ng(k) - 1;
    qgbas(k) = qgend(k-1) + 1;     qgend(k) = qgbas(k) + ng(k) - 1;
  end
end
% contract active energy quantities
pcbas = qgend(nc+1) + 1;
pcend = pcbas + ng(1) - 1;
% upward and downward P reserves
rPpbas = pcend + 1;                rPpend = rPpbas + ng(1) - 1; % positive P reserve
rPmbas = rPpbas + ng(1);           rPmend = rPmbas + ng(1) - 1; % downward P reserve

% if reactive stuff has prices, generate pointers to reactive contract
% quantities, reserves and deviations in the vector of optimization variables
if HAVE_Q
  qcbas = rPmend + 1;            qcend = qcbas + ng(1) - 1;
  rQpbas = qcend + 1;            rQpend = rQpbas + ng(1) - 1; % upward Q reserve
  rQmbas = rQpend + 1;           rQmend = rQmbas + ng(1) - 1; % downward Q reserve
  nvars = rQmend;       % last variable when reactive portion is considered
else
  nvars = rPmend;       % last variable when reactive portion not considered.
end

% starting energy level for each day for short term storage units
if HAVE_DAYS
    s0bas(1) = nvars + 1;       s0end(1) = s0bas(1) + nsts - 1;
    for d = 2:ndays
        s0bas(d) = s0end(d-1) + 1;  s0end(d) = s0bas(d) + nsts - 1;
    end
    nvars = s0end(ndays);
end

% build constraint indexing scheme for nonlinear constraints
pmsmbas(1) = 1;                    pmsmend(1) = pmsmbas(1) + nb(1) - 1;
qmsmbas(1) = pmsmbas(1)+ nb_total; qmsmend(1) = qmsmbas(1) + nb(1) - 1;
sfbas(1) = 2*nb_total + 1;         sfend(1) = sfbas(1) + nl(1) - 1;
stbas(1) = sfbas(1) + nl_total;    stend(1) = stbas(1) + nl(1) - 1;
if nc > 0
  for k = 2:nc+1
    pmsmbas(k) = pmsmend(k-1) + 1; pmsmend(k) = pmsmbas(k) + nb(k) - 1;
    qmsmbas(k) = qmsmend(k-1) + 1; qmsmend(k) = qmsmbas(k) + nb(k) - 1;
    sfbas(k) = sfend(k-1) + 1;     sfend(k) = sfbas(k) + nl(k) - 1;
    stbas(k) = stend(k-1) + 1;     stend(k) = stbas(k) + nl(k) - 1;
  end
end

% SECTION 4: LINEAR CONSTRAINTS FOR ACTIVE POWER RESERVES AND INCREMENTS

aug_gen_stat = [ ones(ng(1),1) gen_stat]; % table includes base case column
                                          % as opposed to gen_stat

% Start defining the constraints.  First define ng(1) upward P reserve variables
% rpp <= rppmax, but do it backwards so we get multipliers with
% correct sign; MINOS writes a Lagrangian F(X) - lambda^T * G(X),
% so we want the lower limit to be usually binding to get positive
% multipliers, so we rewrite this as -rppmax <= -rpp
A1 = sparse((1:ng(1))', (rPpbas:rPpend)', -ones(ng(1),1), ng(1), nvars);
%l1 = -min(PositiveActiveReserveQuantity, gen(:, RAMP_10)) / baseMVA;
l1 = -PositiveActiveReserveQuantity / baseMVA;
u1 = zeros(ng(1),1);
lc1bas = stend(nc+1) + 1;               % linear constraint set 1 start index
lc1end = lc1bas + ng(1)-1;              % linear constraint set 1 end index
% Now define ng(1) downward P reserve variables via rpm <= rpmmax,
% or, for the reasons explained above, as -rpmmax <= -rpm
A2 = sparse((1:ng(1))', (rPmbas:rPmend)', -ones(ng(1),1), ng(1), nvars);
%l2 = -min(NegativeActiveReserveQuantity, gen(:, RAMP_10)) / baseMVA;
l2 = -NegativeActiveReserveQuantity / baseMVA;
u2 = zeros(ng(1), 1);
lc2bas = lc1end + 1;
lc2end = lc2bas + ng(1) - 1;
% alpha-controlled equality of Pi0 and Pc
if FORCE_PC_EQ_P0
  alpha = 0;
else
  alpha = 2*max(max(abs(gen(:,[PMIN PMAX QMIN QMAX]))))/baseMVA+10;
end
A4 = sparse([1:ng(1) 1:ng(1)]', [pgbas(1):pgend(1) pcbas:pcend]', ...
            [ones(ng(1),1); -ones(ng(1),1)],   ng(1), nvars);
l4 = -alpha*ones(ng(1),1);
u4 = alpha*ones(ng(1),1);
lc4bas = lc2end + 1;
lc4end = lc4bas + ng(1) - 1;
% bounds on Pc
A4B = sparse([1:ng(1)]', [pcbas:pcend]', ones(ng(1),1), ng(1), nvars);
l4B = ActiveContractMin / baseMVA;
u4B = ActiveContractMax / baseMVA;
lc4Bbas = lc4end + 1;
lc4Bend = lc4Bbas + ng(1) - 1;
% capacity limits
if ~isempty(caplim_map)
  nclb = size(caplim_map, 1);     %% number of caplim bounds
  A4C = sparse(nclb, nvars);
  A4C(:, (rPpbas:rPpend)') = caplim_map;
%   [ii, jj, ss] = find(caplim_map);
%   A4C = sparse(ii, jj + rPpbas-1, ss, nclb, nvars);
  if isempty(caplim_min)
    l4C = -Inf(nclb, 1);
  else
    l4C = caplim_min / baseMVA;
  end
  if isempty(caplim_max)
    u4C = Inf(nclb, 1);
  else
    u4C = caplim_max / baseMVA;
  end
else
  nclb = 0;
  A4C = sparse(nclb, nvars);
  l4C = [];
  u4C = [];
end
lc4Cbas = lc4Bend + 1;
lc4Cend = lc4Cbas + nclb - 1;

% The Pg minus positive reserve variables (times availability factors) are
% smaller than Pc in all flows. Pg - avail_fac * rPp - Pc <= 0.
% Note that these are the constraints that set
% the shadow price on the upward reserve variables
A6 = sparse(0,0); l6 = []; u6 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A6 = [ A6 ;
      sparse( [1:ng(k)               1:ng(k)           1:ng(k)]', ...
              [(pgbas(k):pgend(k))'; rPpbas-1+ii;      pcbas(1)-1+ii ], ...
              [ones(ng(k),1);       -avail_fac(ii, k); -ones(ng(k),1) ], ...
              ng(k), nvars) ];
  l6 = [ l6 ; -Inf(ng(k), 1) ];
  u6 = [ u6 ; zeros(ng(k), 1) ];
end
lc6bas = lc4Cend + 1;
lc6end = lc6bas + size(A6,1) - 1;
% dispatches are greater than percentage of Gmax (proxy for an aggregate Pmin)
% 0 <= Pg - beta(Pc + avail_fac * rPp)
A7B = sparse(0,0); l7B = []; u7B = [];
for k = 1:nc+1
  genk = Augmented_gen(sum(ng(1:k-1)) + (1:ng(k)), :);
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  % only for gens with PMIN >= 0, no loads, no neg part of storage
  jj = find(genk(:, PMIN) >= 0);
  idx = [pgbas(k):pgend(k)]';
  ngk = length(jj);
  A7B = [ A7B ;
      sparse( [  1:ngk      1:ngk           1:ngk]', ...
            [ idx(jj);      pcbas-1+ii(jj); rPpbas-1+ii(jj) ], ...
            [ ones(ngk,1);  -PminFactor(ii(jj));  -PminFactor(ii(jj)) .* avail_fac(ii(jj), k) ], ...
             ngk, nvars) ];
  l7B = [ l7B ;
         zeros(ngk,1) ];
  u7B = [ u7B;
         Inf(ngk, 1) ];
end
lc7Bbas = lc6end + 1;
lc7Bend = lc7Bbas + size(A7B,1) - 1;

% The Pg plus negative reserve variables (times availability factors) are
% larger than Pc in all flows. Pg + avail_fac * rPm - Pc >= 0.
% Note that these are the constraints that set
% the shadow prices on the downward reserve variables.
A9 = sparse(0,0); l9 = []; u9 = [];
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
  A9 = [ A9 ;
      sparse( [1:ng(k)               1:ng(k)           1:ng(k)]', ...
              [(pgbas(k):pgend(k))'; rPmbas-1+ii;      pcbas(1)-1+ii ], ...
              [ones(ng(k),1);        avail_fac(ii, k); -ones(ng(k),1) ], ...
              ng(k), nvars) ];
  l9 = [ l9 ; zeros(ng(k), 1) ];
  u9 = [ u9 ; Inf(ng(k), 1) ];
end
lc9bas = lc7Bend + 1;
lc9end = lc9bas + size(A9,1) - 1;
% (deleted) the difference between the injection and the contract
% is equal to the inc minus the dec: Pik - Pci = dPp - dPm
A10 = sparse(0,0); l10 = []; u10 = [];
lc10bas = lc9end + 1;
lc10end = lc10bas + size(A10,1) - 1;

% SECTION 5: LINEAR CONSTRAINTS FOR REACTIVE POWER RESERVES AND INCREMENTS

if HAVE_Q
  % Start defining the constraints.  First define ng(1) upward Q reserve variables
  % rqp <= rqpmax, but do it backwards so we get multipliers with
  % correct sign; MINOS writes a Lagrangian F(X) - lambda^T * G(X),
  % so we want the lower limit to be usually binding to get positive
  % multipliers, so we rewrite this as -rqpmax <= -rqp
  A11 = sparse((1:ng(1))', (rQpbas:rQpend)', -ones(ng(1),1), ng(1), nvars);
  l11 = -PositiveReactiveReserveQuantity / baseMVA;
  u11 = zeros(ng(1),1);
  lc11bas = lc10end + 1;        % linear constraint set 11 start index
  lc11end = lc11bas + ng(1)-1;  % linear constraint set 11 end index
  % Now define ng(1) downward Q reserve variables via rqm <= rqmmax,
  % or, for the reasons explained above, as -rqmmax <= -rqm.
  A12 = sparse((1:ng(1))', (rQmbas:rQmend)', -ones(ng(1),1), ng(1), nvars);
  l12 = -NegativeReactiveReserveQuantity / baseMVA;
  u12 = zeros(ng(1), 1);
  lc12bas = lc11end + 1;
  lc12end = lc12bas + ng(1) - 1;
  % alpha-controlled equality of Qi0 and Qc; alpha was computed earlier
  A14 = sparse([1:ng(1) 1:ng(1)]', [qgbas(1):qgend(1) qcbas:qcend]', ...
              [ones(ng(1),1); -ones(ng(1),1)],   ng(1), nvars);
  l14 = -alpha*ones(ng(1),1);
  u14 = alpha*ones(ng(1),1);
  lc14bas = lc12end + 1;
  lc14end = lc14bas + ng(1) - 1;
  % bounds on Qc
  A14B = sparse([1:ng(1)]', [qcbas:qcend]', ones(ng(1),1), ng(1), nvars);
  l14B = ReactiveContractMin / baseMVA;
  u14B = ReactiveContractMax / baseMVA;
  lc14Bbas = lc14end + 1;
  lc14Bend = lc14Bbas + ng(1) - 1;
  % (deleted) incs non-negative
  A15 = sparse(0,0); l15 = []; u15 = [];
  lc15bas = lc14Bend + 1;
  lc15end = lc15bas + size(A15, 1) - 1;
  % The Qg minus positive reserve variables (times availability factors) are
  % smaller than Qc in all flows. Qg - avail_fac * rQp - Qc <= 0.
  % Note that these are the constraints that set
  % the shadow price on the upward reserve variables
  A16 = sparse(0,0); l16 = []; u16 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A16 = [ A16 ;
        sparse( [1:ng(k)               1:ng(k)           1:ng(k)]', ...
                [(qgbas(k):pgend(k))'; rQpbas-1+ii;      qcbas(1)-1+ii ], ...
                [ones(ng(k),1);       -avail_fac(ii, k); -ones(ng(k),1) ], ...
                ng(k), nvars) ];
    l16 = [ l16 ; -Inf(ng(k), 1) ];
    u16 = [ u16 ; zeros(ng(k), 1) ];
  end
  lc16bas = lc15end + 1;
  lc16end = lc16bas + size(A16,1) - 1;
  % (deleted) decs non-negative
  A18 = sparse(0,0); l18 = []; u18 = [];
  lc18bas = lc16end + 1;
  lc18end = lc18bas + size(A18, 1) - 1;
  % The Qg plus negative reserve variables (times availability factors) are
  % larger than Qc in all flows. Qg + avail_fac * rQm - Qc >= 0.
  % Note that these are the constraints that set
  % the shadow prices on the downward reserve variables.
  A19 = sparse(0,0); l19 = []; u19 = [];
  for k = 1:nc+1
    ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
    A19 = [ A19 ;
        sparse( [1:ng(k)               1:ng(k)           1:ng(k)]', ...
                [(qgbas(k):pgend(k))'; rQmbas-1+ii;      qcbas(1)-1+ii ], ...
                [ones(ng(k),1);        avail_fac(ii, k); -ones(ng(k),1) ], ...
                ng(k), nvars) ];
    l19 = [ l19 ; zeros(ng(k), 1) ];
    u19 = [ u19 ; Inf(ng(k), 1) ];
  end
  lc19bas = lc18end + 1;
  lc19end = lc19bas + size(A19,1) - 1;
  % (deleted) the difference between the injection and the contract
  % is equal to the inc minus the dec: Qik - Qci = dQp - dQm
  A20 = sparse(0,0); l20 = []; u20 = [];
  lc20bas = lc19end + 1;
  lc20end = lc20bas + size(A20,1) - 1;
else
  A11 = sparse(0,0); A12 = A11; A14 = A11; A14B = A11;
  A15 = A11; A16 = A11; A18 = A11; A19 = A11; A20 = A11;
  l11 = []; l12 = []; l14 = []; l14B = []; l15 = [];
  l16 = []; l18 = []; l19 = []; l20 = [];
  u11 = []; u12 = []; u14 = []; u14B = []; u15 = [];
  u16 = []; u18 = []; u19 = []; u20 = [];
  lc11bas = lc10end + 1;   lc11end = lc11bas - 1;
  lc12bas = lc11bas;       lc12end = lc11end;
  lc14bas = lc11bas;       lc14end = lc11end;
  lc14Bbas = lc11bas;      lc14Bend = lc11end;
  lc15bas = lc11bas;       lc15end = lc11end;
  lc16bas = lc11bas;       lc16end = lc11end;
  lc18bas = lc11bas;       lc18end = lc11end;
  lc19bas = lc11bas;       lc19end = lc11end;
  lc20bas = lc11bas;       lc20end = lc11end;
end

% total output constraints
if ~isempty(toc_map)
  A21 = sparse(0,0); l21 = []; u21 = [];
  ntoc = length(toc_max);
  for i = 1:ntoc
    for k = 1:nc+1
      ii = find(aug_gen_stat(:, k) > 0);    % which gens active in kth flow
      if k == 1
        p = prob0;
        A21tmp = sparse(1, nvars);
      else
        jj = find(contab(:, CT_LABEL) == clist(k-1));
        p = contab(jj(1), CT_PROB);
      end
      if size(toc_map, 3) == 1      % single map for all k
        k3 = 1;
      else                          % map changes as fcn of k
        k3 = k;
      end
      kk = find(toc_map(i, ii, k3));    % which of these are in constraint i, scenario k
      jj = (pgbas(k):pgend(k))';        % active gen dispatches in kth flow
      A21tmp = A21tmp + ...
          sparse( 1, jj(kk), p * toc_map(i, ii(kk), k3)' .* toc_coeff(ii(kk), toc_type(i)), 1, nvars);
    end
    A21 = [ A21 ; A21tmp ];
  end
  l21 = toc_min / baseMVA;
  u21 = toc_max / baseMVA;
else
  A21 = sparse(0,0); l21 = []; u21 = [];
end
lc21bas = lc20end + 1;
lc21end = lc21bas + size(A21,1) - 1;

if HAVE_DAYS
    % "reserve" capacity equal for postive & negative portion of storage units
    % rpp(p) - rpm(m) = 0
    irpp = (rPpbas:rPpend)';
    irpm = (rPmbas:rPmend)';
    A22 = sparse([1:nsts 1:nsts]', [irpp(istsp); irpm(istsm)], ...
                [ones(nsts,1); -ones(nsts,1)], nsts, nvars);
    l22 = zeros(nsts, 1);
    u22 = l22;
    
    k = 0;
    A23 = sparse(0, nvars); l23 = []; u23 = [];
    ARp = sparse((1:nsts)', irpp(istsp), tsts, nsts, nvars);
    for d = 1:ndays
        hrs = days{d, 1};   % contab indices for each hour of day d
        % fraction of year represented by day d
        d_frac = 0;
        for h = 1:length(hrs)
            if d == 1 && h == 1
                tmp = prob0;
            else
                tmp = contab(hrs(h) == contab(:, CT_LABEL), CT_PROB);
            end
            d_frac = d_frac + tmp(1);
        end

        % state of charge bounds on s(h), state of charge at end of hour h
        % s(h) = s(h-1) - hh(d,h) * (eff * d(d,h) + p(d,h))
        % s(h) = s0 - sum(hh(d,j) * (eff * d(d,j) + p(d,j)))
        %             j=1..h
        As0 = sparse((1:nsts)', s0bas(d):s0end(d), ones(nsts, 1), nsts, nvars);
        Ads = sparse(nsts, nvars);  %% coefficients for change in s through hour h
        for h = 1:length(hrs)
            k = k + 1;
            ii = find(aug_gen_stat(:, k) > 0); % which gens active in kth flow
            g2p = zeros(ng(1), 1);
            g2p(ii) = (1:ng(1))';   % map from gen index to index into pgbas:pgend
            if d == 1 && h == 1
                prob = prob0;
            else
                jj = find(contab(:, CT_LABEL) == hrs(h));
                prob = contab(jj(1), CT_PROB);
            end

            Hdh = days{d, 2} * prob / d_frac;   % real-life hours rep by day d, hour h

            %% change in hour h, -Hdh * (eff * d + p)
            Ads = Ads + Hdh * ...
                sparse([1:nsts 1:nsts]', pgbas(k)-1+g2p([istsp;istsm]), ...
                    [ones(nsts, 1); xi], nsts, nvars);
            
            % state of charge lower bound on s(h) : 0 <= s(h)
            A23 = [A23; As0-Ads];
            l23 = [l23; zeros(nsts, 1)];
            u23 = [u23; Inf(nsts, 1)];
            
            % state of charge upper bound on s(h) : s(h) - tsts * Rp <= 0
            A23 = [A23; As0 - Ads - ARp];
            l23 = [l23; -Inf(nsts, 1)];
            u23 = [u23; zeros(nsts, 1)];
        end
        
        % sum of charge + discharge over day = 0
        A23 = [A23; Ads];
        l23 = [l23; zeros(nsts, 1)];
        u23 = [u23; zeros(nsts, 1)];
    end
else
    A22 = sparse(0,0); l22 = []; u22 = [];
    A23 = sparse(0,0); l23 = []; u23 = [];
end
lc22bas = lc21end + 1;
lc22end = lc22bas + size(A22,1) - 1;
lc23bas = lc22end + 1;
lc23end = lc23bas + size(A23,1) - 1;

% I know, stacking sparse matrices row-wise...
A = [ A1; A2; A4; A4B; A4C; A6; A7B; A9; A10; A11; A12; A14; A14B; A15; A16; A18; A19; A20; A21; A22; A23];
l = [ l1; l2; l4; l4B; l4C; l6; l7B; l9; l10; l11; l12; l14; l14B; l15; l16; l18; l19; l20; l21; l22; l23];
u = [ u1; u2; u4; u4B; u4C; u6; u7B; u9; u10; u11; u12; u14; u14B; u15; u16; u18; u19; u20; u21; u22; u23];

% SECTION 6: Form generalized cost

% Now add cost on reserves; make the cost THE initial linear combination
% in MATPOWER's general cost formulation. This adds a lot of zeros to
% N, but it is easy... change to sparse N add-up of coefficients
% later for efficiency...
Cr = zeros(nvars, 1);
Cr(rPpbas:rPpend) = PositiveActiveReservePrice;
Cr(rPmbas:rPmend) = NegativeActiveReservePrice;
if HAVE_Q
  Cr(rQpbas:rQpend) = PositiveReactiveReservePrice;
  Cr(rQmbas:rQmend) = NegativeReactiveReservePrice;
end
if HAVE_QUADRATIC       % we have quadratic terms
  h = zeros(nvars, 1);
  h(rPpbas:rPpend) = PositiveActiveReservePrice2;
  h(rPmbas:rPmend) = NegativeActiveReservePrice2;
  if HAVE_Q
    h(rQpbas:rQpend) = PositiveReactiveReservePrice2;
    h(rQmbas:rQmend) = NegativeReactiveReservePrice2;
  end
end
if HAVE_QUADRATIC       % we have quadratic terms
  % let N be a sparse identity just to keep it easy to index things
  if isempty(N)
    N = speye(nvars, nvars);
    fparm = ones(nvars, 1) * [1 0 0 1];
    Cw = baseMVA * Cr;
    H = sparse(1:nvars, 1:nvars, 2 * baseMVA^2 * h, nvars, nvars);
  else
    N = [ N;  speye(nvars, nvars) ];
    fparm = [ fparm;  ones(nvars, 1) * [1 0 0 1] ];
    Cw = [ Cw;  baseMVA * Cr ];
    H = [ H sparse(size(H,1),nvars);
          sparse(nvars,size(H,2)) sparse(1:nvars, 1:nvars, 2 * baseMVA^2 * h, nvars, nvars) ];
  end
else                        % no quadratic terms
  % put the linear terms in a single row in N
  if isempty(N)
    N = baseMVA * sparse(Cr');    % make it a sparse row
    fparm = [ 1 0 0 1];           % linear identity fn, no dead zone
    Cw = 1;                       % times 1
    H  = sparse(1,1);             % no quadratic term
  else
    N = [ N;  baseMVA * sparse(Cr') ];
    fparm = [ fparm;  1 0 0 1 ];
    Cw = [ Cw;  1 ];
    H = [ H  sparse(size(H,1),1);  sparse(1,size(H,2)+1 ) ];
  end
end

% SECTION 7: Call solver
mpc = struct(...
    'baseMVA',  baseMVA, ...
    'bus',      Augmented_bus, ...
    'gen',      Augmented_gen, ...
    'branch',   Augmented_branch, ...
    'gencost',  Augmented_gencost, ...
    'A',        A, ...
    'l',        l, ...
    'u',        u, ...
    'N',        N, ...
    'fparm',    fparm, ...
    'H',        H, ...
    'Cw',       Cw ...
);
if ~isempty(iflims)
  mpc.if.map  = Augmented_ifmap;
  mpc.if.lims = Augmented_iflims;
  mpc = toggle_iflims(mpc, 'on');
end
if ~isempty(softlims)
  mpc.softlims.idx = Augmented_softlimsmap;
  mpc.softlims.cost = Augmented_softlimscost;
  %mpc = rmfield(mpc, 'softlims');
  mpc = toggle_softlims(mpc, 'on');
end
if ~isempty(dcline)
  mpc.dcline = Augmented_dcline;
  mpc = toggle_dcline(mpc, 'on');
end
[r, success] = opf(mpc, mpopt);

[buso, geno, brancho, f, info, et] = ...
    deal(r.bus, r.gen, r.branch, r.f, r.raw.info, r.et);
if isfield(r.raw, 'dg')
    [g, jac] = deal(r.raw.g, r.raw.dg);
else
    g = [];
    jac = [];
end
if dc
  pimul = [ r.lin.mu.l.Pmis - r.lin.mu.u.Pmis;
            zeros(size(r.bus, 1), 1);
            -r.branch(:, MU_SF) * baseMVA;
            -r.branch(:, MU_ST) * baseMVA;
            r.lin.mu.l.usr  - r.lin.mu.u.usr    ];
else
  pimul = [ r.mu.nln.l      - r.mu.nln.u;
            r.lin.mu.l.usr  - r.lin.mu.u.usr ];
end

if success && (OUT_ALL > 0)
  printpf(r, 1, mpopt);
end

% SECTION 8: Preliminary unpacking of results

% create an xr that doesn't have any y variables in it
xr = [  buso(:, VA)*pi/180;
        buso(:, VM);
        geno(:, PG)/baseMVA;
        geno(:, QG)/baseMVA;
        r.var.val.z
    ];

% Pick out reserve variables from optimization vector (can be
% larger than actual reserve needed if reserve cost was zero).
% Thus result could be different from Gmax - Pc or Pc - Gmin;
% so it must be verified (see Gmin, Gmax, Qmin, Qmax below).
% For non-unity availability factors, R does not equal
% Gmax - Pc, so we have to use the optimization variables directly
PReservePlus = baseMVA * xr(rPpbas:rPpend);
PReserveMinus = baseMVA * xr(rPmbas:rPmend);
if HAVE_Q
  QReservePlus = baseMVA * xr(rQpbas:rQpend);
  QReserveMinus = baseMVA * xr(rQmbas:rQmend);
else
  QReservePlus = [];
  QReserveMinus = [];
end

% Pick out dispatches, reserve, and some multipliers, in internal order
Pik = zeros(ng(1), nc+1);
Qik = zeros(ng(1), nc+1);
lamPik = zeros(ng(1), nc+1);
lamQik = zeros(ng(1), nc+1);
Pik(:, 1) = geno(1:ng(1), PG);  % the base case dispatch
Qik(:, 1) = geno(1:ng(1), QG);
Pc = baseMVA * xr(pcbas:pcend); % the contract active dispatch
if HAVE_Q
  Qc = baseMVA * xr(qcbas:qcend); % the contracted reactive dispatch
end
lamPik(:, 1) = buso(geno(1:ng(1), GEN_BUS), LAM_P); % base case multipliers
lamQik(:, 1) = buso(geno(1:ng(1), GEN_BUS), LAM_Q);
lamRPp_GT_dPpik = zeros(ng(1), nc+1);
lamRPp_GT_dPpik(:, 1) = -pimul(lc6bas-1+(1:ng(1))) / baseMVA;
lamRPm_GT_dPmik = zeros(ng(1), nc+1);
lamRPm_GT_dPmik(:, 1) = pimul(lc9bas-1+(1:ng(1))) / baseMVA;

lamRPplusmax = max(0,pimul(lc1bas:lc1end)/baseMVA); % upward reserve -ramp <= -Rp+ <= 0
%lamRPplusmax = pimul(lc1bas:lc1end)/baseMVA; % upward reserve -ramp <= -Rp+ <= 0
%lamRPplusmin = min(0,pimul(lc1bas:lc1end)/baseMVA); % negative means 0 limit binding
lamRPminusmax = max(0,pimul(lc2bas:lc2end)/baseMVA); % downward reserve -ramp <= -Rp- <= 0
%lamRPminusmax = pimul(lc2bas:lc2end)/baseMVA; % downward reserve -ramp <= -Rp- <= 0
%lamRPminusmin = min(0,pimul(lc2bas:lc2end)/baseMVA);
lam_alpha = -pimul(lc4bas-1+(1:ng(1))) / baseMVA;
lam_Pc = -pimul(lc4Bbas-1+(1:ng(1))) / baseMVA;

if HAVE_Q
  lamRQp_GT_dQpik = zeros(ng(1), nc);
  lamRQp_GT_dQpik(:, 1) = -pimul(lc16bas-1+(1:ng(1))) / baseMVA;
  lamRQm_GT_dQmik = zeros(ng(1), nc);
  lamRQm_GT_dQmik(:, 1) = pimul(lc19bas-1+(1:ng(1))) / baseMVA;
  lamQikplus = zeros(ng(1), nc+1);
  lamRQplusmax = pimul(lc11bas:lc11end) / baseMVA;
  lamRQminusmax = pimul(lc12bas:lc12end) / baseMVA;
  lam_alphaQ = -pimul(lc14bas-1+(1:ng(1))) / baseMVA;
  lam_Qc = -pimul(lc14Bbas-1+(1:ng(1))) / baseMVA;
else
  lamRQp_GT_dQpik = [];
  lamRQm_GT_dQmik = [];
  lamQikplus = [];
  lamQikminus = [];
  lamRQplusmax = [];
  lamRQminusmax = [];
  lam_alphaQ = [];
end

% Pick out dispatches and several multipliers for each contingency.
for k = 1:nc+1
  ii = find(aug_gen_stat(:, k) > 0);
  Pik(ii, k) = geno(sum(ng(1:k-1))+1:sum(ng(1:k)), PG); % this is how to pick from augm Gen table, works for k=1
  Qik(ii, k) = geno(sum(ng(1:k-1))+1:sum(ng(1:k)), QG); % because sum(empty) = 0.
  lamPik(ii, k) = buso(geno(sum(ng(1:k-1))+1:sum(ng(1:k)), GEN_BUS), LAM_P);
  lamQik(ii, k) = buso(geno(sum(ng(1:k-1))+1:sum(ng(1:k)), GEN_BUS), LAM_Q);  
  lamRPp_GT_dPpik(ii, k) = -pimul((lc6bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  lamRPm_GT_dPmik(ii, k) = pimul((lc9bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  if HAVE_Q
    lamRQp_GT_dQpik(ii, k) = -pimul((lc16bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
    lamRQm_GT_dQmik(ii, k) = pimul((lc19bas-1+sum(ng(1:k-1))) + (1:ng(k))) / baseMVA;
  end
end

% Find Gmax/Gmin, but do it only over periods where the generator
% was explicitly committed.  Thus can't do it with max()
% and min() over the matrices directly.
Gmin = zeros(ng(1), 1);
Gmax = zeros(ng(1), 1);
for i=1:ng(1)
  ii = find([ 1 gen_stat(i, :)] > 0);
  Gmin(i) = min(Pik(i, ii));
  Gmax(i) = max(Pik(i, ii));
end
Qmin = zeros(ng(1), 1);   % these are interesting even when Q is not
Qmax = zeros(ng(1), 1);   % considered in the market/cost
for i=1:ng(1)
  ii = find([ 1 gen_stat(i, :)] > 0);
  Qmin(i) = min(Qik(i, ii));
  Qmax(i) = max(Qik(i, ii));
end

% short term storage
if nsts
    % consolidate, then remove rows for negative gen part of short term storage
    PReservePlus(istsm) = [];                       % delete negative gen part
    PReserveMinus(istsp) = PReserveMinus(istsm);    % copy negative gen part
    PReserveMinus(istsm) = [];                      % delete negative gen part
    
    Pik(istsp, :) = Pik(istsp, :) + Pik(istsm, :);
    Pik(istsm, :) = [];                             % delete negative gen part
    Qik(istsp, :) = Qik(istsp, :) + Qik(istsm, :);
    Qik(istsm, :) = [];                             % delete negative gen part

    Gmin(istsp) = Gmin(istsm);                      % copy negative gen part
    Gmin(istsm) = [];                               % delete negative gen part
    Gmax(istsm) = [];                               % delete negative gen part
    Qmin(istsm) = [];                               % delete negative gen part
    Qmax(istsm) = [];                               % delete negative gen part

    lamPik(istsm, :) = [];                          % delete negative gen part
    lamQik(istsm, :) = [];                          % delete negative gen part

    Pc(istsm) = [];                                 % delete negative gen part

    lamRPp_GT_dPpik(istsp, :) = lamRPp_GT_dPpik(istsp, :) + lamRPp_GT_dPpik(istsm, :);
    lamRPp_GT_dPpik(istsm, :) = [];                 % delete negative gen part
    lamRPm_GT_dPmik(istsp, :) = lamRPm_GT_dPmik(istsp, :) + lamRPm_GT_dPmik(istsm, :);
    lamRPm_GT_dPmik(istsm, :) = [];                 % delete negative gen part

    lamRPplusmax(istsp, :) = lamRPplusmax(istsp, :) + lamRPplusmax(istsm, :);
    lamRPplusmax(istsm, :) = [];                    % delete negative gen part
    lamRPminusmax(istsp, :) = lamRPminusmax(istsp, :) + lamRPminusmax(istsm, :);
    lamRPminusmax(istsm, :) = [];                   % delete negative gen part

    lam_alpha(istsp, :) = lam_alpha(istsp, :) + lam_alpha(istsm, :);
    lam_alpha(istsm, :) = [];                       % delete negative gen part

    lam_Pc(istsp, :) = abs(lam_Pc(istsp, :)) + abs(lam_Pc(istsm, :));
    lam_Pc(istsm, :) = [];                          % delete negative gen part

    % extract s0 for each day
    s0 = zeros(nsts, ndays);
    for d = 1:ndays
        s0(:, d) = xr(s0bas(d):s0end(d)) * baseMVA;
    end
% s0
end

% Sum some multipliers over flows
sumlamP = sum(lamPik, 2);
sumlamQ = sum(lamQik, 2);
sumlamRPp = sum(lamRPp_GT_dPpik, 2);
sumlamRPm = sum(lamRPm_GT_dPmik, 2);
if HAVE_Q
  sumlamRQp = sum(lamRQp_GT_dQpik, 2);
  sumlamRQm = sum(lamRQm_GT_dQmik, 2);
end

% lamP0_div_prob0 = lamPik(:, 1) / prob0;


% SECTION 9: Create output results

results = struct;

% base flow results
results.base.bus = buso(1:nb(1), :);
branchwork = brancho(1:nl(1), :);
genwork = geno(1:ng(1), :);
genwork(:, VG) = buso(geno(1:ng(1), GEN_BUS), VM);
[results.base.bus, genwork, branchwork] = ...
   int2ext(i2e, results.base.bus, genwork, branchwork);
results.base.branch = branch_original;
tmp = zeros(length(base_off_branch), length([PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
results.base.branch(base_off_branch, [PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = tmp;
results.base.branch(base_on_branch, :) = branchwork;
results.base.gen = gen_original;
tmp = zeros(length(base_off_gen), length([PG;QG;MU_PMAX;MU_PMIN;MU_QMAX;MU_QMIN]));
results.base.gen(base_off_gen, [PG;QG;MU_PMAX;MU_PMIN;MU_QMAX;MU_QMIN]) = tmp;
if nsts
    genwork(istsp, PG) = genwork(istsp, PG) + genwork(istsm, PG);
    genwork(istsp, MU_PMAX) = genwork(istsp, MU_PMAX) + genwork(istsm, MU_PMIN);
    genwork(istsp, MU_PMIN) = genwork(istsp, MU_PMIN) + genwork(istsm, MU_PMAX);
    genwork(istsm, :) = [];                         % delete negative gen part
end
results.base.gen(base_on_gen, :) = genwork;
if ~isempty(dcline)
  dclineo = r.dcline;
  dclinework = dclineo(1:ndc(1), :);
  dclinework(:, c.F_BUS) = i2e( dclinework(:, c.F_BUS) );
  dclinework(:, c.T_BUS) = i2e( dclinework(:, c.T_BUS) );
  results.base.dcline = dcline_original;
  tmp = zeros(length(base_off_dcline), length([c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
  results.base.dcline(base_off_dcline, [c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = tmp;
  results.base.dcline(base_on_dcline, :) = dclinework;
end
if ~isempty(iflims)
  results.base.if = iflims;
  results.base.if.P = r.if.P(1:nif(1));
  results.base.if.mu.l = r.if.mu.l(1:nif(1));
  results.base.if.mu.u = r.if.mu.u(1:nif(1));
end
if ~isempty(softlims)
    results.base.softlims = softlims; % If user-specified fiels are added
    results.base.softlims.idx = softlims.idx;
    results.base.softlims.cost = softlims.cost(:, :, 1);
    results.base.softlims.overload = r.softlims.overload(Augmented_softlimsmap(1:nsoftlims(1)));
    results.base.softlims.ovl_cost = r.softlims.ovl_cost(Augmented_softlimsmap(1:nsoftlims(1)));
    %results.base.softlims.mu.u = r.softlims.mu.u(1:nsoftlims(1));
end

% post-contingency flows results
kk = 1;  % index over actual considered contingencies
for k=1:length(clist_original) % index over original list of contingencies, 
  clabel = clist_original(k);  % including those thrown out over uncommit issues
  if any(clabel == clist) % Was this contingency actually included in analysis?
    % pull kth contingency bus data from augmented flow
    buswork = buso((1:nb(kk+1))+sum(nb(1:kk)), :);
    buswork(:, BUS_I) = (1:nb(kk+1))';
    % branch back-transformation in two steps: first to coincide with
    % "branch", which already does not have uncommitted branches on input, then
    % to coincide with branch_original, which is the original input data
    results.cont(k).branch = branch_original;
    branchwork = branch;
    branchwork1 = brancho((1:nl(kk+1))+sum(nl(1:kk)), :); %  pull from augmented
    branchwork1(:, [F_BUS;T_BUS]) = ...
        branchwork1(:, [F_BUS;T_BUS]) - sum(nb(1:kk));    % renumber FROM, TO bus #
    branchwork(branch_stat(:, kk) ~= 0, :) = branchwork1; % 1st back-transform step;
    ii = find(~branch_stat(:, kk));  % now zero out fields in inactive branches
    branchwork(ii, [BR_STATUS;PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = ...
      zeros(length(ii), length([BR_STATUS;PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
    % still have to do 2nd back-transformation but since we must first
    % reorder gens, we turn our attention to them...
    results.cont(k).gen = gen_original;
    genwork = gen;
    ii = find(gen_stat(:, kk) > 0);
    genwork(ii, :) = geno((1:ng(kk+1))+sum(ng(1:kk)), :); % pull from augm flow
    genwork(ii, GEN_BUS) = genwork(ii, GEN_BUS) - sum(nb(1:kk)); % backmap bus #
    genwork(ii, VG) = buso(geno((1:ng(kk+1))+sum(ng(1:kk))-1, GEN_BUS), VM);
    ii = find(gen_stat(:, kk) <= 0); % zero out fields for gen taken out in contingency..
    genwork(ii, [GEN_STATUS;PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]) = ...
      zeros(length(ii), length([GEN_STATUS;PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]));
    if nsts
        genwork(istsp, PG) = genwork(istsp, PG) + genwork(istsm, PG);
        genwork(istsp, MU_PMAX) = genwork(istsp, MU_PMAX) + genwork(istsm, MU_PMIN);
        genwork(istsp, MU_PMIN) = genwork(istsp, MU_PMIN) + genwork(istsm, MU_PMAX);
        genwork(istsm, :) = [];                     % delete negative gen part
    end
    % renumber back the buses
    [buswork, genwork, branchwork] = ...
        int2ext(i2e, buswork, genwork, branchwork);
    % and stick into results struct, zeroing out appropriate fields for
    % gens and branches that were uncommitted on input...
    results.cont(k).bus = buswork;
    results.cont(k).branch(base_on_branch, :) = branchwork;
    tmp = zeros(length(base_off_branch), length([PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]));
    results.cont(k).branch(base_off_branch, [PF;QF;PT;QT;MU_SF;MU_ST;MU_ANGMIN;MU_ANGMAX]) = tmp;
    results.cont(k).gen(base_on_gen, :) = genwork;
    tmp = zeros(length(base_off_gen), length([PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]));
    results.cont(k).gen(base_off_gen, [PG;QG;MU_PMIN;MU_PMAX;MU_QMIN;MU_QMAX]) = tmp;
    if ~isempty(dcline)
      % dcline back-transformation in two steps: first to coincide with
      % "dcline", which already does not have uncommitted dclines on input, then
      % to coincide with dcline_original, which is the original input data
      results.cont(k).dcline = dcline_original;
      dclinework = dcline;
      dclinework1 = dclineo((1:ndc(kk+1))+sum(ndc(1:kk)), :); %  pull from augmented
      dclinework1(:, [c.F_BUS;c.T_BUS]) = ...
          dclinework1(:, [c.F_BUS;c.T_BUS]) - sum(nb(1:kk));    % renumber FROM, TO bus #
%       dclinework(dcline_stat(:, kk) ~= 0, :) = dclinework1; % 1st back-transform step;
%       ii = find(~dcline_stat(:, kk));  % now zero out fields in inactive dclines
%       dclinework(ii, [c.BR_STATUS;c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = ...
%         zeros(length(ii), length([c.BR_STATUS;c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
      dclinework = dclinework1;
      % still have to do 2nd back-transformation but since we must first
      % reorder gens, we turn our attention to them...
      % renumber back the buses
      dclinework(:, c.F_BUS) = i2e( dclinework(:, c.F_BUS) );
      dclinework(:, c.T_BUS) = i2e( dclinework(:, c.T_BUS) );
      % and stick into results struct, zeroing out appropriate fields for
      % dclines that were uncommitted on input...
      results.cont(k).dcline(base_on_dcline, :) = dclinework;
      tmp = zeros(length(base_off_dcline), length([c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]));
      results.cont(k).dcline(base_off_dcline, [c.PF;c.PT;c.QF;c.QT;c.VF;c.VT;c.MU_PMIN;c.MU_PMAX;c.MU_QMINF;c.MU_QMAXF;c.MU_QMINT;c.MU_QMAXT]) = tmp;
    end
    if ~isempty(iflims)
      results.cont(k).if = iflims;
      results.cont(k).if.P = r.if.P((1:nif(kk+1))+sum(nif(1:kk)));
      results.cont(k).if.mu.l = r.if.mu.l((1:nif(kk+1))+sum(nif(1:kk)));
      results.cont(k).if.mu.u = r.if.mu.u((1:nif(kk+1))+sum(nif(1:kk)));
    end
    if ~isempty(softlims)
      results.cont(k).softlims = softlims;
      results.cont(k).softlims.idx = softlims.idx;
      results.cont(k).softlims.cost = softlims.cost(:, :, k+1);
      results.cont(k).softlims.overload = r.softlims.overload(Augmented_softlimsmap((1:nsoftlims) + (nsoftlims * k)));
      results.cont(k).softlims.ovl_cost = r.softlims.ovl_cost(Augmented_softlimsmap((1:nsoftlims) + (nsoftlims * k)));
      %results.cont(k).if.mu.u = r.if.mu.u((1:nif(kk+1))+sum(nif(1:kk)));
    end
    kk = kk + 1;
  else % this contingency was deleted from the local list as a result of some
    results.cont(k).bus = [];  % equipment being uncommited on input, and
    results.cont(k).branch = []; % this contingency attempting to take it out,
    results.cont(k).gen = [];    % so we just insert empty contingency flow data
  end
end

% Compute sum of bus lambdas over contingencies plus base case
sum_bus_lam_p = buso(1:nb(1), LAM_P);
sum_bus_lam_q = buso(1:nb(1), LAM_Q);
for k = 2:nc+1
  sum_bus_lam_p = sum_bus_lam_p + buso((1:nb(k))+sum(nb(1:k-1)), LAM_P);
  sum_bus_lam_q = sum_bus_lam_q + buso((1:nb(k))+sum(nb(1:k-1)), LAM_Q);
end
results.energy.prc.sum_bus_lam_p = sum_bus_lam_p;
results.energy.prc.sum_bus_lam_q = sum_bus_lam_q;

ng_original = size(gen_original, 1);

results.energy.sum_muPmax = results.base.gen(:, MU_PMAX);
results.energy.sum_muPmin = results.base.gen(:, MU_PMIN);
for k=1:nc
  if ~isempty(results.cont(k).gen)
    results.energy.sum_muPmax = results.energy.sum_muPmax + ...
                                results.cont(k).gen(:, MU_PMAX);
    results.energy.sum_muPmin = results.energy.sum_muPmin + ...
                                results.cont(k).gen(:, MU_PMIN);
  end
end

% The following output data is all in ng_original x (nc+1) vectors; nc
% is actual number of contingencies considered in this run
tmp = zeros(ng_original, nc+1);
results.reserve.mu.Rp_pos = tmp; % upward reserve for active power
results.reserve.mu.Rp_neg = tmp; % downward reserve for active power
results.reserve.mu.Rp_pos(base_on_gen, :) = lamRPp_GT_dPpik;
results.reserve.mu.Rp_neg(base_on_gen, :) = lamRPm_GT_dPmik;
if HAVE_Q
  results.reserve.mu.Rq_pos = tmp;
  results.reserve.mu.Rq_neg = tmp;
  results.reserve.mu.Rq_pos(base_on_gen, :) = lamRQp_GT_dQpik;
  results.reserve.mu.Rq_neg(base_on_gen, :) = lamRQm_GT_dQmik;
end

results.energy.Pc = zeros(ng_original, 1);
results.energy.Pc(base_on_gen) = Pc;
if HAVE_Q
  results.energy.Qc = zeros(ng_original, 1);
  results.energy.Qc(base_on_gen) = Qc;
end
tmp = zeros(ng_original, 1);
results.energy.Gmax = tmp;
results.energy.Gmin = tmp;
results.energy.Qmax = tmp;
results.energy.Qmin = tmp;
results.energy.Gmax(base_on_gen) = Gmax;
results.energy.Gmin(base_on_gen) = Gmin;
results.energy.Qmax(base_on_gen) = Qmax;
results.energy.Qmin(base_on_gen) = Qmin;

% The following output data is in (ng_original x 1) vectors.
tmp = zeros(ng_original, 1);
results.energy.mu.alphaP = tmp;
results.energy.mu.Pc = tmp;
results.reserve.qty.Rp_pos = tmp;
results.reserve.qty.Rp_neg = tmp;
results.reserve.prc.Rp_pos = tmp;
results.reserve.prc.Rp_neg = tmp;
results.reserve.mu.Rpmax_pos = tmp;
results.reserve.mu.Rpmax_neg = tmp;
results.energy.mu.alphaP(base_on_gen) = lam_alpha;
results.energy.mu.Pc(base_on_gen) = lam_Pc;
results.reserve.qty.Rp_pos(base_on_gen) = PReservePlus;
results.reserve.qty.Rp_neg(base_on_gen) = PReserveMinus;
results.reserve.prc.Rp_pos(base_on_gen) = sumlamRPp;
results.reserve.prc.Rp_neg(base_on_gen) = sumlamRPm;
results.reserve.mu.Rpmax_pos(base_on_gen) = lamRPplusmax;
results.reserve.mu.Rpmax_neg(base_on_gen) = lamRPminusmax;
if HAVE_Q
  results.energy.mu.alphaQ = tmp;
  results.energy.mu.Qc = tmp;
  results.reserve.qty.Rq_pos = tmp;
  results.reserve.qty.Rq_neg = tmp;
  results.reserve.prc.Rq_pos = tmp;
  results.reserve.prc.Rq_neg = tmp;
  results.reserve.mu.Rqmax_pos = tmp;
  results.reserve.mu.Rqmax_neg = tmp;
  results.energy.mu.alphaQ(base_on_gen) = lam_alphaQ;
  results.energy.mu.Qc(base_on_gen) = lam_Qc;
  results.reserve.qty.Rq_pos(base_on_gen) = QReservePlus;
  results.reserve.qty.Rq_neg(base_on_gen) = QReserveMinus;
  results.reserve.prc.Rq_pos(base_on_gen) = sumlamRQp;
  results.reserve.prc.Rq_neg(base_on_gen) = sumlamRQm;
  results.reserve.mu.Rqmax_pos(base_on_gen) = lamRQplusmax;
  results.reserve.mu.Rqmax_neg(base_on_gen) = lamRQminusmax;
end
if ~isempty(caplim_map)
    if nsts
        caplim_map(:, istsm) = [];                  % delete negative gen part
    end
    results.caplim.map = caplim_map;
    results.caplim.min = caplim_min;
    results.caplim.max = caplim_max;
    results.caplim.qty = A4C * xr * baseMVA;
    results.caplim.mu  = -pimul(lc4Cbas:lc4Cend) / baseMVA;
end
if ~isempty(toc_map)
    if nsts
        toc_map(:, istsm, :) = [];                  % delete negative gen part
        toc_coeff(istsm, :) = [];                   % delete negative gen part
    end
    results.total_output.map = toc_map;
    results.total_output.min = toc_min;
    results.total_output.max = toc_max;
    results.total_output.coeff = toc_coeff;
    results.total_output.type = toc_type;
    results.total_output.qty = A21 * xr * baseMVA;
    results.total_output.mu = -pimul(lc21bas:lc21end) / baseMVA;
end
if HAVE_DAYS
    results.days = days;
end
if nsts
    results.short_term_storage.s0 = zeros(nsts0, ndays);
    results.short_term_storage.s0(base_sts_on, :) = s0;
end

% raw OPF results
results.success = success;
results.opf_results = r;

% Eval marginal cost of all injections
if HAVE_Q
  xx = [ geno(:, PG); geno(:, QG)];
else
  xx = geno(:, PG);
end
df_dPg = margcost(Augmented_gencost, xx);
if HAVE_Q
  df_dQg = df_dPg(ng_total+1:2*ng_total);
  df_dPg = df_dPg(1:ng_total);
end
%df_dPg_out = zeros(ng_original, 1);
%df_dPg_out(base_on_original) = df_dPg;  NO

% Eval marginal cost of offers at Pc quantities
if HAVE_Q
  xx = [ Pc; Qc];
else
  xx = Pc;
end
if nsts
    gen(istsm, :) = [];                             % delete negative gen part
    gencost(istsm, :) = [];                         % delete negative gen part
end
df_dPcQc = margcost(gencost, xx);
df_dPc = zeros(ng_original, 1);
df_dPc(base_on_gen) = df_dPcQc;
if HAVE_Q
  df_dQc = zeros(ng_original, 1);
  df_dQc(base_on_gen) = df_dPcQc(ng(1)+1:2*ng(1));
end

tmp = zeros(ng_original, 1);
sumlamP_out = tmp;
sumlamP_out(base_on_gen) = sumlamP;
sum_genbus_lam_p_out = tmp;
sum_genbus_lam_p_out(base_on_gen) = sum_bus_lam_p(gen(:, GEN_BUS));

if success && OUT_ALL
  ngencol = 6; % number of generators per print table (max # columns)
  fprintf('\n');
  fprintf('================================================================================\n')
  fprintf('  E4ST Results\n');
  fprintf('================================================================================\n')
  fprintf('\n');
  fprintf('Total expected cost: $%g \n', f);
  fprintf('Total R+: %g \n', sum(PReservePlus));
  fprintf('Total R-: %g \n', sum(PReserveMinus));
  fprintf('Generator summary\n');
  fprintf('\n');
  for i=1:ngencol:ng_original
    igp = i:min(i+ngencol-1, ng_original);
    fprintf('           ');
    fprintf('   G%2i   ', igp);
    fprintf('\n');
    fprintf('================================================================================\n')
    fprintf('    Pc   |');
    fprintf(' %7.2f ', results.energy.Pc(igp));
    fprintf('\n');
    fprintf('   Gmin  |');
    fprintf(' %7.2f ', results.energy.Gmin(igp));
    fprintf('\n');
    fprintf('   Gmax  |');
    fprintf(' %7.2f ', results.energy.Gmax(igp));
    fprintf('\n');
%    fprintf('    P0   |');
%    fprintf(' %7.2f ', results.base.gen(igp, PG));
%    fprintf('\n');
    fprintf(' P(i, 0) |');
    fprintf(' %7.2f ', results.base.gen(igp, PG));
    fprintf('\n');
    for k = 1:nc
%      fprintf(' P(i,%2i) |', k);   % this prints k
      fprintf(' P(i,%2i) |', clist(k)); % but this prints the contingency label
      fprintf(' %7.2f ', results.cont(k).gen(igp, PG));
      fprintf('\n');
    end
    fprintf('  C''(Pc) |');
    fprintf(' %7.2f ', df_dPc(igp));
    fprintf('\n');
    fprintf(' Sum Lpg |');
    fprintf(' %7.2f ', sumlamP_out(igp));
    fprintf('\n');
    fprintf(' Lp_bus  |');
    fprintf(' %7.2f ', sum_genbus_lam_p_out(igp));
    fprintf('\n');
    fprintf('sumMuPmax|');
    fprintf(' %7.2f ', results.energy.sum_muPmax(igp));
    fprintf('\n');
    fprintf('sumMuPmin|');
    fprintf(' %7.2f ', results.energy.sum_muPmin(igp));
    fprintf('\n');
    fprintf('   Rp+   |');
    fprintf(' %7.2f ', results.reserve.qty.Rp_pos(igp));
    fprintf('\n');
    fprintf('  $Rp+   |');
    fprintf(' %7.2f ', results.reserve.prc.Rp_pos(igp));
    fprintf('\n');
    fprintf('muRp+box |');
    fprintf(' %7.2f ', results.reserve.mu.Rpmax_pos(igp));
    fprintf('\n');
    fprintf('   Rp-   |');
    fprintf(' %7.2f ', results.reserve.qty.Rp_neg(igp));
    fprintf('\n');
    fprintf('  $Rp-   |');
    fprintf(' %7.2f ', results.reserve.prc.Rp_neg(igp));
    fprintf('\n');
    fprintf('muRp-box |');
    fprintf(' %7.2f ', results.reserve.mu.Rpmax_neg(igp));
    fprintf('\n');
  end
end


%sumlamP
%lam_alpha
%lamP0_div_prob0
%sumlamQ
%sumlamRPp
%sumlamRPp_div_prob0 = sumlamRPp / prob0
%sumlamRPm
%sumlamRPm_div_prob0 = sumlamRPm / prob0
%lamRPplusmax
%lamRPminusmax
%sum_bus_lam_p
%sum_bus_lam_q


%disp('Gradient of cost wrt all injections')
%df_dPg

%disp('marginal costs at Pc/Qc quantities:')
%df_dPcQc

