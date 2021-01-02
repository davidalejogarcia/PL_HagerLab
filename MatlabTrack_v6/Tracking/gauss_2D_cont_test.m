%% Test how Gaussian step function looks

x_center = 0;
y_center = 0;

x = (-5:.2:5);
y = (-5:.2:5);


%% Test implementation
parFit(1) = 100;
parFit(2) = 1.5;
parFit(3) = 0;

parMod{1} = 1;
parMod{2} = x;
parMod{3} = y;
parMod{4} = x_center;
parMod{5} = y_center;

I_test = gauss_2D_cont_fun(parFit,parMod);

figure
h=surfc(x,y,I_test);
set(h,'EdgeColor','k');
%set(gca,'XGrid','off')
%set(gca,'YGrid','off')
%set(gca,'ZGrid','off')
axis([-5 5 -5 5 -0 100])

%% Test fit
[parFit ssr I_fit] = gauss_2D_cont_fit(I_test,parMod,parFit);
