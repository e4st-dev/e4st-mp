function mpc = applyUplift(mpc, pol_val, idx_gen)

if ~isfield(mpc, 'uplift')
    mpc.uplift = zeros(length(idx_gen), 1);
end

mpc.uplift(idx_gen, 1) = pol_val;
