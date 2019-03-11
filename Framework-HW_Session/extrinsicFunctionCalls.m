function r = ext_MTG(x_hat, param)
r = zeros(10, 1);
coder.extrinsic('myTargetGenerator');
r = myTargetGenerator( x_hat, param);
end

function x_hat = ext_MSE(u, y, param)
x_hat = zeros(16,1);
coder.extrinsic('myStateEstimator');
x_hat = myStateEstimator(u, y, param);
end

function u = ext_MPC(r, x_hat, param)
u = zeros(2, 1);
coder.extrinsic('myMPController');
u = myMPController(r, x_hat, param);
end
