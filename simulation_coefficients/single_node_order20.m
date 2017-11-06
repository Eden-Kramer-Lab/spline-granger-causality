function b = single_node_order20
% Returns coefficients for univariate AR(20) model
b = zeros(1,1,20);
b(1,1,:) = [-0.2314 %%%%% spectral peak, high freq
    0.1613
    0.1405
    0.0741
    -0.0098
    -0.0836
    -0.1193
    -0.1048
    -0.0585
    0.0011
    0.0559
    0.0874
    0.0837
    0.0614
    0.0289
    -0.0050
    -0.0316
    -0.0423
    -0.0285
    0.0185];


end