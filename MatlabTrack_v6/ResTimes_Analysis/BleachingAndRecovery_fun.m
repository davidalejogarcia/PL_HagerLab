function esp = BleachingAndRecovery_fun(coef,x,val0)

k_bl = coef(1);
k_rec = coef(2);
k_mot = coef(3);
a = coef(4);

M = [-k_bl-k_mot, k_rec, k_mot; k_bl, -k_rec, 0; k_mot, 0, -k_mot];

[V,D] = eig(M);

c = V\val0;
% 
% c1 = (val0(1) - val0(2)*(V(1,2)/V(2,2)))/(V(1,1) - (V(1,2)*V(2,1)/V(2,2)));
% 
% c2 = val0(2)-c1*V(2,1)/V(2,2);

esp = c(1)*V(1,1)*exp(D(1,1)*x) + c(2)*V(1,2)*exp(D(2,2)*x) + c(3)*V(1,3)*exp(D(3,3)*x);
esp = a*esp;

% esp = exp((-k_bl + k_rec)*x);
% esp = a*esp;

