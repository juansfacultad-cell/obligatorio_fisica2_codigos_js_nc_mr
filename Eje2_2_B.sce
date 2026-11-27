// 2.2.b - Bz(0,z) sobre el eje, z en [0, L/2] (+ extension opcional)
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

// ---------- Calculo sobre el eje ----------
Nz = 1000;
z  = linspace(0, L/2, Nz);
Bz_ax = zeros(1, Nz);
for i = 1:Nz
    Bz_ax(i) = (mu0*I*a^2/2) * sum( 1 ./ (a^2 + (z(i) - zk).^2).^(3/2) );
end

scf(3); clf();
plot(z, Bz_ax, "b-");
xlabel("z (m)"); ylabel("B_z(0, z) (T)");
title("Campo sobre el eje del solenoide largo");
legend("B_z numerico", "in_lower_left"); xgrid();
xs2png(3, "fig_2_2_b.png");

printf("Bz(0)   = %.5e T\n", Bz_ax(1));
printf("Bz(L/2) = %.5e T\n", Bz_ax($));

// ---------- Extension opcional mas alla del extremo ----------
// Si la usas, indicalo y justificalo en el informe (ver como decae el campo fuera).
extender = %t;
if extender then
    zext = linspace(0, 1.2*L/2, Nz);
    Bext = zeros(1, Nz);
    for i = 1:Nz
        Bext(i) = (mu0*I*a^2/2) * sum( 1 ./ (a^2 + (zext(i) - zk).^2).^(3/2) );
    end
    scf(4); clf();
    plot(zext, Bext, "b-");
    plot([L/2 L/2], [0 max(Bext)], "k--");
    xlabel("z (m)"); ylabel("B_z(0, z) (T)");
    title("Campo sobre el eje, extendido mas alla del extremo");
    legend("B_z numerico", "extremo z = L/2"); xgrid();
    xs2png(4, "fig_2_2_b_extendido.png");
end
