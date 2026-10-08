function amb = ref_ambiguity(ref, gap)
%REF_AMBIGUITY Smallest Delta E from each row to any row whose label differs by more than GAP.
%   Small values mean the colour alone cannot tell this structure from a different one.
if nargin < 2, gap = ref.label.gap; end
n = numel(ref.value);
amb = nan(n, 1);
for s = 1:1000:n                                  % chunks keep memory at 1000 x n
    e = min(s + 999, n);
    D = sqrt(max(sum(ref.Jab(s:e, :) .^ 2, 2) + sum(ref.Jab .^ 2, 2)' - 2 * ref.Jab(s:e, :) * ref.Jab', 0));
    D(abs(ref.value(s:e) - ref.value') <= gap) = inf;
    m = min(D, [], 2);
    m(~isfinite(m)) = nan;
    amb(s:e) = m;
end
end
