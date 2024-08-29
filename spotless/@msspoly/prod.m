function q = prod(p,d)
    if nargin < 2,
        d = 1;
        if size(p,d) == 1
            d = 3-d;
        end
    end
    
    switch d,
      case 1,
      case 2,
        p = p';
      otherwise,
        error('d must be in {1,2}.');
    end
    
    q = ones(1,size(p,2));

    for i = 1:size(p,1)
        q = q.*indexinto(p,i,':');
    end

    if d == 2
        q = q';
    end
end