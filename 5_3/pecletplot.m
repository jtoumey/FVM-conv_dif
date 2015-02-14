rho = 1;
u = -2.5;
F_w = rho*u;
phi_W = 1;
phi_P = 2;
gamma = 0.1;

L = 0:.01:1;
Pe_w = rho*u.*L/gamma;




q_w = (0.5*(1 + 2./Pe_w)*phi_W + 0.5*(1 - 2./Pe_w)*phi_P);

plot(L,q_w)