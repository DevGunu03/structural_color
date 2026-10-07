function [idx, dist, altVal, altDist] = nearest_ref(X, ref, valid, gap)
%NEAREST_REF Nearest reference row (J'a'b') for each pixel, exhaustive and chunked.
%   With GAP (MoO3), also the best row whose value differs from the best match by > GAP.
R = ref.Jab;
rr = sum(R .^ 2, 2)';
n = size(X, 1);
idx = zeros(n, 1);  dist = zeros(n, 1);
wantAlt = nargout > 2;
altVal = nan(n, 1);  altDist = nan(n, 1);
chunk = 20000;
for s = 1:chunk:n
    e = min(s + chunk - 1, n);
    x = X(s:e, :);
    d2 = sum(x .^ 2, 2) + rr - 2 * (x * R');
    [m, j] = min(d2, [], 2);
    idx(s:e) = j;
    dist(s:e) = sqrt(max(m, 0));
    if wantAlt
        v = valid(s:e);
        near = abs(ref.value' - ref.value(j)) <= gap;
        d2(near) = inf;
        [m2, j2] = min(d2, [], 2);
        ok = v & isfinite(m2);
        rows = (s:e)';
        altVal(rows(ok)) = ref.value(j2(ok));
        altDist(rows(ok)) = sqrt(max(m2(ok), 0));
    end
end
end
