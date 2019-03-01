function x_hat = ext_MSE(u, y, param)
x_hat = zeros(16,1);
coder.extrinsic('myStateEstimator');
x_hat = myStateEstimator(u, y, param);
end

