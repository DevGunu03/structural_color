function [idx, dist, altVal, altDist, uniq] = nearest_ref(X, ref, valid, gap, tol, sigma)
%NEAREST_REF Nearest reference row (J'a'b') for each pixel, exhaustive and chunked.
%   With GAP, also the best row whose label differs from the best match by more than GAP
%   (interference colours repeat, so a different structure can have almost the same colour).
%   With TOL and SIGMA, also the uniqueness of the match: the share of the weight
%   w_k = prior_k exp(-d_k^2 / 2 sigma^2) that lies within TOL of the best match (or in its class).
%   Same as match_statistics in python/tcd/reliability.py.
R = ref.Jab;
rr = sum(R .^ 2, 2)';
n = size(X, 1);
idx = zeros(n, 1);  dist = zeros(n, 1);
wantAlt = nargout > 2;
wantU = nargout > 4;
altVal = nan(n, 1);  altDist = nan(n, 1);  uniq = nan(n, 1);
if wantU, prior = label_prior(ref.value)'; end
chunk = 20000;
for s = 1:chunk:n
    e = min(s + chunk - 1, n);
    x = X(s:e, :);
    d2 = max(sum(x .^ 2, 2) + rr - 2 * (x * R'), 0);
    [m, j] = min(d2, [], 2);
    idx(s:e) = j;
    dist(s:e) = sqrt(m);
    if wantAlt
        v = valid(s:e);
        best = ref.value(j);
        rows = (s:e)';
        if wantU
            w = prior .* exp(-(d2 - m) / (2 * sigma ^ 2));
            u = sum(w .* same_structure(ref, ref.value', best, tol), 2) ./ sum(w, 2);
            uniq(rows(v)) = u(v);
        end
        d2(abs(ref.value' - best) <= gap) = inf;
        [m2, j2] = min(d2, [], 2);
        ok = v & isfinite(m2);
        altVal(rows(ok)) = ref.value(j2(ok));
        altDist(rows(ok)) = sqrt(m2(ok));
    end
end
end

function w = label_prior(v)
% Candidate weights uniform in the label: half the distance to each neighbour, capped at twice
% the typical step (same as label_prior in reliability.py).
[s, order] = sort(v(:));
gaps = diff(s);
typical = median(gaps(gaps > 0));
if isempty(typical) || isnan(typical), typical = 1; end
left = min([typical; gaps], 2 * typical);
right = min([gaps; typical], 2 * typical);
w = zeros(size(s));
w(order) = max((left + right) / 2, 1e-12 * typical);
w = w / sum(w);
end
