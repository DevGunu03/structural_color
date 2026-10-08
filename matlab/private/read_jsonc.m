function def = read_jsonc(file)
%READ_JSONC Decode a JSON file that may contain // comments (outside strings).
text = fileread(file);
out = blanks(numel(text));
n = 0;
i = 1;
inStr = false;
while i <= numel(text)
    c = text(i);
    if inStr
        n = n + 1; out(n) = c;
        if c == '\' && i < numel(text)
            i = i + 1;
            n = n + 1; out(n) = text(i);
        elseif c == '"'
            inStr = false;
        end
    elseif c == '"'
        inStr = true;
        n = n + 1; out(n) = c;
    elseif c == '/' && i < numel(text) && text(i + 1) == '/'
        while i <= numel(text) && text(i) ~= newline
            i = i + 1;
        end
        continue
    else
        n = n + 1; out(n) = c;
    end
    i = i + 1;
end
try
    def = jsondecode(out(1:n));
catch err
    error('tcd:system', '%s: not valid JSON after removing // comments: %s', file, err.message);
end
end
