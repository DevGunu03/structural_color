function M = fit_correction(measured, target, mode, ridge)
%FIT_CORRECTION 3x3 matrix with target ~ measured * M' (linear RGB, one row per anchor).
%   'diagonal' : per-channel gains; 'matrix' : full matrix ridge-regularised toward the gains;
%   'auto'     : diagonal for one anchor, matrix otherwise.
if nargin < 4, ridge = 1e-2; end
D = diag(sum(measured .* target, 1) ./ sum(measured .^ 2, 1));
if strcmp(mode, 'diagonal') || (strcmp(mode, 'auto') && size(measured, 1) == 1)
    M = D;
    return
end
G = measured' * measured;
lam = ridge * trace(G) / 3;
M = ((G + lam * eye(3)) \ (measured' * target + lam * D'))';
end
