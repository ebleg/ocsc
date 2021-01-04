function out = vaf(x, xest)
    out = max(0, 100*(1 - var(x - xest)/var(x)));
end