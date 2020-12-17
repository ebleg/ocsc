function [re, im] = gen_nyquist(H)
    assert(all(size(H) == [2 2]), 'H must be 2x2.');
    
    arg = eye(2) + H; % This function will only work for 2x2
    det_arg = minreal(arg(1,1)*arg(2,2) - arg(2,1)*arg(1,2), [], false);
    [re, im] = nyquist(det_arg); % Use the determinant
    re = squeeze(re); im = squeeze(im);
end

