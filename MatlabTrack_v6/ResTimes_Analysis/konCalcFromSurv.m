function [k_on, tau_search] = konCalcFromSurv(SurvFitPar)


%Calculates the pseudo-on rate (k_on*Seq) from an SMT experiment using the
%2-component exponential fit of the survival histogram.

%non-specific-bound fraction of all bound mols
Fns = SurvFitPar(1,5);
Fns_err = SurvFitPar(2,5);
%specific-bound fraction of all bound mols
Fs = 1 - Fns;
Fs_err = SurvFitPar(2,5);
%Bound-Fraction
BF = SurvFitPar(1,6);
BF_err = SurvFitPar(2,6);

%specific binding
tau_s = 1/SurvFitPar(1,4);
tau_s_err = (tau_s/SurvFitPar(1,4))*SurvFitPar(2,4);

%non-specific binding
tau_ns = 1/SurvFitPar(1,3);
tau_ns_err = (tau_ns/SurvFitPar(1,3))*SurvFitPar(2,3);

%average binding
tau_bar = Fs*tau_s + Fns*tau_ns;
tau_bar_err1 = Fs*tau_s*sqrt((Fs_err/Fs).^2 + (tau_s_err/tau_s).^2);
tau_bar_err2 = Fns*tau_ns*sqrt((Fns_err/Fns).^2 + (tau_ns_err/tau_ns).^2);

tau_bar_err = sqrt(tau_bar_err1.^2 + tau_bar_err2.^2);

%k_on
k_on = BF./((1+Fns/Fs)*(1-BF)*tau_bar);
% k_on_den1 = (1+Fns/Fs);
% k_on_den2 = (1-BF)*tau_bar;

k_on_err_den1 = (1 + Fns/Fs)*sqrt((Fns_err/Fns).^2 + (Fs_err/Fs).^2);
k_on_err_den2 = ((1-BF)*tau_bar)*sqrt((BF_err/BF).^2 + (tau_bar_err/tau_bar).^2);
% k_on_err_den = ((1+Fns/Fs)*(1-BF)*tau_bar)*sqrt(
k_on_err = k_on*sqrt((BF_err/BF).^2 + (k_on_err_den1/(1+Fns/Fs)).^2 + (k_on_err_den2/((1-BF)*tau_bar)).^2);

k_on = [k_on, k_on_err];

% Search Time

tau3D = tau_bar*(1-BF)/BF;

tau3D_err = tau3D*sqrt((tau_bar_err/tau_bar).^2 + (BF_err/(1-BF)).^2 + (BF_err/BF).^2);

Ntrials = 1/Fs;
Ntrials_err = Ntrials*(Fs_err/Fs);

tau_search = Ntrials*tau3D + (Ntrials - 1)*tau_ns;
tau_search_err1 = Ntrials*tau3D*sqrt((Ntrials_err/Ntrials).^2 + (tau3D_err/tau3D).^2);
tau_search_err2 = (Ntrials - 1)*tau_ns*sqrt((Ntrials_err/(Ntrials - 1)).^2 + (tau_ns_err/tau_ns).^2);
tau_search_err = sqrt(tau_search_err1.^2 + tau_search_err2.^2);

tau_search = [tau_search, tau_search_err];

