function tf = same_structure(ref, a, b, tol)
%SAME_STRUCTURE True where labels A and B count as the same structure (same class, or within TOL).
if ref.label.classes
    tf = floor(a + 0.5) == floor(b + 0.5);
else
    tf = abs(a - b) <= tol;
end
end
