function n_eff = maxwell_garnett(n_host, n_incl, f)
%MAXWELL_GARNETT Effective index of spheres (volume fraction F) in a host medium.
if f < 0 || f > 1
    error('tcd:fill', 'Fill fraction %g outside [0, 1]', f);
end
eh = n_host .^ 2;
ei = n_incl .^ 2;
n_eff = sqrt(eh .* (2 * (1 - f) * eh + (1 + 2 * f) * ei) ./ ((2 + f) * eh + (1 - f) * ei));
end
