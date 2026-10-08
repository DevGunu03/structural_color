function n = ema_mix(rule, n_host, n_incl, f)
%EMA_MIX Effective index of a two-phase mixture (inclusion volume fraction F in a host).
%   'maxwell-garnett' : isolated spheres in a host (as in TransferMatrix_packing.m)
%   'bruggeman'       : both phases on an equal footing; of the two roots of
%                       2 e^2 - b e - e_i e_h = 0, b = (3f - 1) e_i + (2 - 3f) e_h, the one with the
%                       larger imaginary part (or, for lossless media, the positive one)
%   'linear'          : volume-averaged permittivity, e = f e_i + (1 - f) e_h
%   Same formulas as python/tcd/materials.py.
if f < 0 || f > 1
    error('tcd:fill', 'Fill fraction %g outside [0, 1]', f);
end
switch rule
    case 'maxwell-garnett'
        n = maxwell_garnett(n_host, n_incl, f);
    case 'bruggeman'
        eh = complex(n_host .^ 2);
        ei = complex(n_incl .^ 2);
        b = (3 * f - 1) * ei + (2 - 3 * f) * eh;
        root = sqrt(b .^ 2 + 8 * ei .* eh);
        e1 = (b + root) / 4;
        e2 = (b - root) / 4;
        pick1 = imag(e1) > imag(e2) + 1e-12 | (abs(imag(e1) - imag(e2)) <= 1e-12 & real(e1) >= real(e2));
        e = e2;
        e(pick1) = e1(pick1);
        n = sqrt(e);
    case 'linear'
        n = sqrt(complex(f * n_incl .^ 2 + (1 - f) * n_host .^ 2));
    otherwise
        error('tcd:ema', 'Unknown mixing rule ''%s'' (maxwell-garnett, bruggeman or linear)', rule);
end
end
