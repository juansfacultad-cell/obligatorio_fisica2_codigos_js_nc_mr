// 2.2.a - Br(r,z) y Bz(r,z) para z = a, 0 < r < 3a
clear; clc; clf();

// ---------- Parametros del solenoide largo ----------
mu0 = 4*%pi*1e-7;
I   = 1;
a   = 0.03;
L   = 100*a;          // 3 m
N   = 20000;
n   = N/L;            // 6666.67 espiras/m
M   = 2048;           // subdivisiones angulares

dz  = L/N;
k   = (0:N-1)';
zk  = -L/2 + (k + 1/2)*dz;     // centros de las espiras (N x 1)

dth   = 2*%pi/M;
theta = ((0:M-1) + 1/2)*dth;   // 1 x M

// ---------- Calculo ----------
z0   = a;
Npts = 90;                              // puntos en r (bajar a 30 para probar rapido)
r    = ((1:Npts) - 0.5)*(3*a/Npts);     // 0 < r < 3a, evita r = a exacto
Br   = zeros(1, Npts);
Bz   = zeros(1, Npts);

blk = 500;                              // espiras por bloque (bajar si falta memoria)
cth = cos(theta);

tic();
for i = 1:Npts
    ri = r(i);
    sBr = 0; sBz = 0;
    for ini = 1:blk:N
        fin = min(ini+blk-1, N);
        K   = fin - ini + 1;
        DZ  = (z0 - zk(ini:fin)) * ones(1, M);
        CT  = ones(K,1) * cth;
        D   = (ri^2 + DZ.^2 + a^2 - 2*a*ri*CT).^(3/2);
        sBr = sBr + sum( a*DZ.*CT ./ D );
        sBz = sBz + sum( (a^2 - a*ri*CT) ./ D );
    end
    Br(i) = mu0*I*dth/(4*%pi) * sBr;
    Bz(i) = mu0*I*dth/(4*%pi) * sBz;
end
printf("Tiempo: %.1f s\n", toc());

// ---------- Graficas ----------
scf(1); clf();
plot(r, Br, "b-");
xlabel("r (m)"); ylabel("B_r (T)");
title("Componente radial B_r(r, z = a)");
legend("B_r numerico", "in_upper_left"); xgrid();
xs2png(1, "fig_2_2_a_Br.png");

scf(2); clf();
plot(r, Bz, "r-");
xlabel("r (m)"); ylabel("B_z (T)");
title("Componente axial B_z(r, z = a)");
legend("B_z numerico", "in_lower_left"); xgrid();
xs2png(2, "fig_2_2_a_Bz.png");

printf("Br(r->0) = %.3e T (esperado ~0 por simetria)\n", Br(1));
printf("Bz(r->0) = %.5e T\n", Bz(1));
