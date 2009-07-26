function [w,u,p,r] = data_processing(thefile, varargin)
% arguments:
% - file to process
% - what to get: real or imag part or both (say 'auto')
% - inferior limit
% - resolution
% - superior limit


optargin = size(varargin, 2);
stdargin = nargin - optargin;

data = load(thefile, '-ascii');

binf = 1;
res  = 1;
bsup = size(data,1);

if optargin > 0
    % real part, anyone?
    if varargin{1} == 'real'
        w = data(:,1);
        u = data(:,3);
        p = data(:,5);
        r = data(:,7);
    elseif varargin{1} == 'imag'
        w = data(:,2);
        u = data(:,4);
        p = data(:,6);
        r = data(:,8);
    else
        w = data(:,1) + i.*data(:,2);
        u = data(:,3) + i.*data(:,4);
        p = data(:,5) + i.*data(:,6);
        r = data(:,7) + i.*data(:,8);
    end

    if ~isempty(varargin{2})
        binf = varargin{2};
    end

    if ~isempty(varargin{3})
        res = varargin{3};
    end

    if ~isempty(varargin{4})
        bsup = varargin{4};
    end
else
    % on demand…
    for k = 1:optargin
        % TODO
        disp 'ERROR: not implemented yet!'
    end
end

w = w(binf:res:bsup);
u = u(binf:res:bsup);
p = p(binf:res:bsup);
r = r(binf:res:bsup);

