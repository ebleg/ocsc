function yhat = simtheta(theta, u, y, Asize, Bsize, Csize, Dsize, Ksize, xsize)
    [A, B, C, D, K, x0] = theta2matrices(theta, Asize, Bsize, ...
                                              Csize, Dsize, Ksize, xsize);
                                          
    yhat = simsystem(A, [B K], C, [D zeros(l, l)], x0, [u y]);
end