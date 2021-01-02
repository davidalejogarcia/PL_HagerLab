function [parFit ssr I_fit] = gauss_2D_cont_fit(I_exp,parMod,parStart)

fun_dum = @(parFit) sum ( sum( ( gauss_2D_cont_fun(parFit,parMod)-I_exp).^2 ) );


parFit = fminsearch(fun_dum,parStart);
I_fit  = gauss_2D_cont_fun(parFit,parMod);
ssr    = sum(sum((I_fit-I_exp).^2));

