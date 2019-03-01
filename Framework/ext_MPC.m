function u = ext_MPC(r, x_hat, param)
u = zeros(2, 1);
coder.extrinsic('myMPController');
u = myMPController(r, x_hat, param);
end
