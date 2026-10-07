function out = mode_filter(labels, sz, nClasses)
%MODE_FILTER Majority filter for an integer label map; -1 (unassigned) is its own class.
[H, W] = size(labels);
r = floor(sz / 2);
ir = [r:-1:1, 1:H, H:-1:H - r + 1];  ir = min(max(ir, 1), H);
ic = [r:-1:1, 1:W, W:-1:W - r + 1];  ic = min(max(ic, 1), W);
k = ones(sz) / sz ^ 2;
best = -inf(H, W);
out = -ones(H, W);
for c = -1:nClasses - 1
    cnt = conv2((labels(ir, ic) == c) * 1.0, k, 'valid');
    take = cnt > best;
    best(take) = cnt(take);
    out(take) = c;
end
end
