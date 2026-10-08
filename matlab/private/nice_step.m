function w = nice_step(span, n)
%NICE_STEP A round bin width (1, 2, 2.5 or 5 x 10^k) giving about N bins over SPAN.
if nargin < 2, n = 12; end
raw = max(span, 1e-12) / n;
p = 10 ^ floor(log10(raw));
m = [1 2 2.5 5 10] * p;
w = m(find(m >= raw - 1e-12, 1));
end
