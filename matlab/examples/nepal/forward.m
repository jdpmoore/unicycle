% This function is called by the copula with sampled model parameters. The 
% physical model is called from here. Here as an example, we call afterslip only
% model. We provided the bounds in config.dat file in logspace, so we convert 
% them to real numbers before passing to physical model.   
% 
% Author : Sagar Masuti.
% -----------------------------------------------------------------------------

function [dsim]=forward(minput)
     minput = 10.^minput;
     dsim=postseismic_copula(minput);
end

