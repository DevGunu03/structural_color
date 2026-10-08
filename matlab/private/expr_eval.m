function v = expr_eval(expr, vars)
%EXPR_EVAL Value of a system-file expression such as 'bead_nm * sqrt(2/3)'.
%   v = expr_eval(expr, vars): EXPR is a number or a char/string; VARS a struct of numbers or
%   equal-length column vectors (constants and sweep parameters). The result broadcasts.
%
%   Same language as python/tcd/expr.py, parsed here (never passed to eval):
%     numbers 2, 0.5, .5, 1e-3; names in VARS and pi; + - * / ^ (or **); unary -; ( );
%     comparisons < <= > >= == != (1 or 0); functions sqrt exp log log10 sin cos tan asin acos
%     atan abs floor ceil round (halves up) and two-argument min(a, b), max(a, b).
%   Precedence, low to high: comparison, + -, * /, unary -, ^ (right-associative).
if isnumeric(expr) || islogical(expr)
    v = double(expr);
    return
end
if isstring(expr), expr = char(expr); end
if ~ischar(expr)
    error('tcd:expr', 'expected a number or an expression, got a %s', class(expr));
end
if nargin < 2, vars = struct(); end
toks = tokenize(expr);
if isempty(toks), error('tcd:expr', 'empty expression'); end
pos = 1;
v = comparison();
if pos <= numel(toks)
    error('tcd:expr', 'unexpected ''%s'' in ''%s''', toks{pos}{2}, expr);
end
if ~isreal(v)
    error('tcd:expr', '''%s'' gives a complex number (square root or log of a negative value?)', expr);
end

    % ---- recursive descent, sharing toks / pos ------------------------------------------
    function t = peek()
        if pos <= numel(toks), t = toks{pos}; else, t = {'end', ''}; end
    end
    function t = take(expected)
        t = peek();
        if nargin > 0 && ~(strcmp(t{1}, 'op') && strcmp(t{2}, expected))
            got = t{2}; if isempty(got), got = 'end of expression'; end
            error('tcd:expr', 'expected ''%s'' but found ''%s'' in ''%s''', expected, got, expr);
        end
        pos = pos + 1;
    end
    function tf = isop(varargin)
        t = peek();
        tf = strcmp(t{1}, 'op') && any(strcmp(t{2}, varargin));
    end
    function a = comparison()
        a = additive();
        if isop('<', '<=', '>', '>=', '==', '!=')
            op = take();
            b = additive();
            switch op{2}
                case '<',  a = double(a < b);
                case '<=', a = double(a <= b);
                case '>',  a = double(a > b);
                case '>=', a = double(a >= b);
                case '==', a = double(a == b);
                otherwise, a = double(a ~= b);
            end
        end
    end
    function a = additive()
        a = term();
        while isop('+', '-')
            op = take();
            b = term();
            if op{2} == '+', a = a + b; else, a = a - b; end
        end
    end
    function a = term()
        a = unary();
        while isop('*', '/')
            op = take();
            b = unary();
            if op{2} == '*', a = a .* b; else, a = a ./ b; end
        end
    end
    function a = unary()
        if isop('-', '+')
            op = take();
            a = unary();
            if op{2} == '-', a = -a; end
        else
            a = power_();
        end
    end
    function a = power_()
        a = atom();
        if isop('^')
            take();
            a = a .^ unary();
        end
    end
    function a = atom()
        t = take();
        switch t{1}
            case 'num'
                a = str2double(t{2});
            case 'name'
                name = t{2};
                if isop('(') && is_function(name)
                    take('(');
                    args = {comparison()};
                    while isop(',')
                        take();
                        args{end + 1} = comparison(); %#ok<AGROW>
                    end
                    take(')');
                    a = call(name, args);
                elseif strcmp(name, 'pi') && ~isfield(vars, 'pi')
                    a = pi;
                elseif isfield(vars, name)
                    a = double(vars.(name));
                else
                    known = strjoin(sort(fieldnames(vars))', ', ');
                    if isempty(known), known = 'none'; end
                    error('tcd:expr', 'unknown name ''%s'' in ''%s'' (known names: %s)', name, expr, known);
                end
            otherwise
                if strcmp(t{2}, '(')
                    a = comparison();
                    take(')');
                else
                    got = t{2}; if isempty(got), got = 'end of expression'; end
                    error('tcd:expr', 'unexpected ''%s'' in ''%s''', got, expr);
                end
        end
    end
    function a = call(name, args)
        need = 1;
        if any(strcmp(name, {'min', 'max'})), need = 2; end
        if numel(args) ~= need
            error('tcd:expr', '%s() takes %d argument(s), got %d in ''%s''', name, need, numel(args), expr);
        end
        x = args{1};
        switch name
            case 'round', a = floor(x + 0.5);           % halves up, as in the Python code
            case 'min',   a = min(x, args{2});
            case 'max',   a = max(x, args{2});
            otherwise,    a = feval(name, x);           % sqrt exp log log10 sin cos tan asin acos atan abs floor ceil
        end
    end
end

function tf = is_function(name)
tf = any(strcmp(name, {'sqrt', 'exp', 'log', 'log10', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan', ...
    'abs', 'floor', 'ceil', 'round', 'min', 'max'}));
end

function toks = tokenize(s)
% Cell array of {kind, text}: kind = 'num' | 'name' | 'op'.
pat = '\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)|([A-Za-z_]\w*)|(\*\*|<=|>=|==|!=|[-+*/^(),<>]))';
toks = {};
s = strtrim(s);
i = 1;
while i <= numel(s)
    e = regexp(s(i:end), ['^' pat], 'end', 'once');
    if isempty(e)
        error('tcd:expr', 'cannot read ''%s'' in expression ''%s''', strtrim(s(i:end)), s);
    end
    txt = strtrim(s(i:i + e - 1));
    if ~isempty(regexp(txt, '^(\d|\.\d)', 'once'))
        toks{end + 1} = {'num', txt}; %#ok<AGROW>
    elseif ~isempty(regexp(txt, '^[A-Za-z_]', 'once'))
        toks{end + 1} = {'name', txt}; %#ok<AGROW>
    else
        if strcmp(txt, '**'), txt = '^'; end
        toks{end + 1} = {'op', txt}; %#ok<AGROW>
    end
    i = i + e;
end
end
