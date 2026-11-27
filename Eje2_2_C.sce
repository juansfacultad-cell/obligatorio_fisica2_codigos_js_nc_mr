// 2.2.c - Comparacion: numerico vs analitico finito vs ideal
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

Nz = 1000;
z  = linspace(0, L/2, Nz);

// Numerico (Biot-Savart, suma sobre espiras)
Bz_ax = zeros(1, Nz);
for i = 1:Nz
    Bz_ax(i) = (mu0*I*a^2/2) * sum( 1 ./ (a^2 + (z(i) - zk).^2).^(3/2) );
end

// Analitico finito
Bfin = (mu0*n*I/2) * ( (z + L/2)./sqrt(a^2 + (z + L/2).^2) ...
                      - (z - L/2)./sqrt(a^2 + (z - L/2).^2) );
// Ideal
Bid = mu0*n*I;

scf(5); clf();
plot(z, Bz_ax, "b-");
plot(z, Bfin, "r--");
plot(z, Bid*ones(1,Nz), "k:");
xlabel("z (m)"); ylabel("B_z(0, z) (T)");
title("Comparacion: numerico, analitico finito e ideal");
legend("Numerico (Biot-Savart)", "Analitico finito", "Ideal mu0 n I", "in_lower_left");
xgrid();
xs2png(5, "fig_2_2_c.png");

// Error relativo respecto al ideal
err = abs(Bz_ax - Bid) / abs(Bid) * 100;

scf(6); clf();
plot(z/L, err, "m-");
plot([0 0.5], [1 1], "k--");
xlabel("z / L"); ylabel("Error relativo (%)");
title("Error relativo del modelo numerico respecto del ideal");
legend("Error relativo", "Umbral 1 %", "in_upper_left"); xgrid();
xs2png(6, "fig_2_2_c_error.png");

// Diferencia numerico vs analitico finito
errNF = abs(Bz_ax - Bfin) ./ abs(Bfin) * 100;
printf("Max diferencia numerico vs analitico finito: %.4f %%\n", max(errNF));

// Tabla de referencia
fr = [0 0.40 0.45 0.49 0.50];
printf("\n   z/L      Bz (T)         Error rel. (%%)\n");
for q = 1:length(fr)
    zq = fr(q)*L;
    Bq = (mu0*I*a^2/2) * sum( 1 ./ (a^2 + (zq - zk).^2).^(3/2) );
    printf("  %.2f   %.5e   %9.3f\n", fr(q), Bq, abs(Bq - Bid)/Bid*100);
end

// Region con error <= 1 %
idx  = find(err <= 1);
zmax = z(idx($));
printf("\nError <= 1 %% hasta z = %.4f m = %.3f L -> region |z| <= %.3f L (%.1f %% de la longitud)\n", ...
        zmax, zmax/L, zmax/L, 2*zmax/L*100);
