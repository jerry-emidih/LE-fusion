function str = literalNum(A, varargin)

    p = inputParser;
    addRequired(p, 'A', @isnumeric)
    addOptional(p, 'pres', 2, @isnumeric)
    parse(p, A, varargin{:})
    
    A = p.Results.A;
    pres = p.Results.pres;
    
    format = sprintf('%%.%dg', pres);

    n = size(A, 1);
    m = size(A, 2);
    if m > 1
        str = join(compose(format, A));
    else
        str = compose(format, A);
    end
    str = replace(str, " ", ", ");
    str = repmat("[", [n, 1]) + str + repmat("]", [n, 1]);
end


