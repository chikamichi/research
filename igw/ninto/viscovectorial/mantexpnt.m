function [expnt, mant] = mantexpnt(arg)
% returns exponant and mantisse
% TODO handle > 1D!

if size(arg,1) == 1 && size(arg,2) == 1
    % scalar
    % nothing to do yet but may be useful someday
elseif size(arg,1) == 1 || size(arg,2) == 1
    % row
    % looking for the max, whatever the sign
    tmin = mantexpnt(min(arg));
    tmax = mantexpnt(max(arg));
    if abs(tmin) > tmax
        arg = min(arg);
    else
        arg = max(arg);
    end
else
    error('does not handle matrix or multidimensionnal arrays yet!');
end

sgn   = sign(arg);
expnt = fix(log10(abs(arg)));
mant  = sgn * 10^(log10(abs(arg))-expnt);

if abs(mant) < 1
    mant = mant * 10;
    expnt = expnt - 1;
end

% if expnt is actually 0
if expnt == -Inf || expnt == Inf
    expnt = 0;
end

