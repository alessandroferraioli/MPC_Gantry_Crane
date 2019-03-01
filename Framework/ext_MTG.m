function r = ext_MTG(x_hat, param)
r = zeros(10, 1);
coder.extrinsic('myTargetGenerator');
r = myTargetGenerator( x_hat, param);
end

