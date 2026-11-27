// ================================================================
// FISICA 2 - Espira circular: apartado 1.1 (a, b y c)
// Campo magnetico por la ley de Biot-Savart.
//
// Convencion: corriente en sentido antihorario vista desde +z.
// Por la regla de la mano derecha, en el centro B apunta hacia +z.
// Unidades SI: a en m, I en A, B en T.
//
// Metodo numerico: regla del punto medio en phi, con M divisiones:
// phi_j = (j + 1/2)*dphi, dphi = 2*pi/M.
// Se agrupan puntos de observacion en bloques para limitar memoria.
// ================================================================

clear;
clc;

// Cerrar figuras anteriores (compatible con distintas versiones de Scilab).
try
    close(winsid());
catch
    try
        delete(winsid());
    catch
    end
end

// ---------------- Parametros fisicos y numericos ----------------
mu0 = 4*%pi*1e-7;     // permeabilidad del vacio (T*m/A)
a   = 0.03;           // radio: 3 cm = 0.03 m
I   = 1.0;             // corriente (A)
Np  = 1000;            // puntos de observacion por recta
blk = 20;              // puntos procesados por bloque

// 1000 puntos estrictamente interiores al intervalo (-5a, 5a).
t = -5*a + (1:Np)'*(10*a/(Np+1));

// ---------------- Campo mediante Biot-Savart -------------------
// Parametrizacion de la espira:
// r'(phi) = (a*cos(phi), a*sin(phi), 0)
// dl' = a*(-sin(phi), cos(phi), 0) dphi
// Las componentes de dl' x (r-r') son:
//   a*z*cos(phi), a*z*sin(phi), a^2-a*(x*cos(phi)+y*sin(phi)).
// Verificacion dimensional: dl' x (r-r') / |r-r'|^3 tiene unidad 1/m;
// mu0*I tiene unidad T*m; por tanto cada contribucion queda en teslas (T).
function [Bx, By, Bz] = campo_espira(xp, yp, zp, a, I, M, blk)
    mu0 = 4*%pi*1e-7;
    n = length(xp);
    dphi = 2*%pi/M;
    phi = ((0:M-1) + 0.5)*dphi;  // puntos medios, vector fila
    cph = cos(phi);
    sph = sin(phi);

    Bx = zeros(n,1);
    By = zeros(n,1);
    Bz = zeros(n,1);
    factor = mu0*I/(4*%pi)*dphi;

    for i = 1:blk:n
        j = min(i+blk-1, n);
        m = j-i+1;

        xi = matrix(xp(i:j), m, 1);
        yi = matrix(yp(i:j), m, 1);
        zi = matrix(zp(i:j), m, 1);
        X = xi*ones(1,M);
        Y = yi*ones(1,M);
        Z = zi*ones(1,M);
        C = ones(m,1)*cph;
        S = ones(m,1)*sph;

        // |r-r'|^3
        R3 = ((X-a*C).^2 + (Y-a*S).^2 + Z.^2).^1.5;

        // Suma numerica (regla del punto medio) de las tres componentes.
        Bx(i:j) = factor*sum(a*Z.*C ./ R3, "c");
        By(i:j) = factor*sum(a*Z.*S ./ R3, "c");
        Bz(i:j) = factor*sum((a^2-a*(X.*C+Y.*S)) ./ R3, "c");
    end
endfunction

// ---------------- Prueba de convergencia -----------------------
// Punto cercano al conductor pero fuera de la espira:
// (x,y,z) = (0.99a, 0, 0). La distancia al hilo es 0.01a.
// Se compara Bz al aumentar la resolucion angular M.
Mconv = [1024, 2048, 4096, 8192];
Bz_anterior = %nan;
mprintf("\nPRUEBA DE CONVERGENCIA en (x,y,z)=(0.99a, 0, 0)\n");
mprintf("a = %.5f m, I = %.3f A\n", a, I);
for k = 1:length(Mconv)
    [bxt, byt, bzt] = campo_espira(0.99*a, 0, 0, a, I, Mconv(k), 1);
    mprintf("M = %6.0f : Bz = %.12e T", Mconv(k), bzt);
    if k > 1 then
        errrel = abs(bzt-Bz_anterior)/abs(bzt);
        mprintf("   diferencia relativa respecto al M anterior = %.3e", errrel);
    end
    mprintf("\n");
    Bz_anterior = bzt;
end

// Comprobacion del sentido y del orden de magnitud en el centro.
// Para la corriente elegida, B(0,0,0) = (mu0*I/(2a))*zhat.
Mcheck = 4096;
[bxc, byc, bzc] = campo_espira(0, 0, 0, a, I, Mcheck, 1);
Bcentro_teorico = mu0*I/(2*a);
mprintf("\nCOMPROBACION EN EL CENTRO DE LA ESPIRA\n");
mprintf("Bz numerico  = %.12e T\n", bzc);
mprintf("Bz teorico   = %.12e T\n", Bcentro_teorico);
mprintf("El signo positivo indica direccion +z para corriente antihoraria vista desde +z.\n\n");

// ---------------- Definicion de las 8 rectas --------------------
// Tipo 1: x fijo, y=0; varia z entre -5a y 5a.
// Tipo 2: y=0, z fijo; varia x entre -5a y 5a.
// ATENCION: las ultimas tres condiciones fijan z, por lo que la
// coordenada libre de esas rectas es x, no z. Por eso sus graficas
// se presentan en funcion de x (m). La consigna parece tener aqui
// una errata cuando dice que las 24 graficas son funcion de z.
//
// Columnas de la matriz: tipo, coordenada fija, divisiones angulares M.
rectas = [ 1, 0,       4096;
           1, a/2,     4096;
           1, 0.99*a, 16384;
           1, 1.01*a, 16384;
           1, 2*a,    4096;
           2, 0.01,   8192;   // z = 0.01 m (no 0.01a)
           2, a,      4096;   // z = 0.03 m
           2, 2*a,    4096 ]; // z = 0.06 m

nombres = ["x=0, y=0";
           "x=a/2, y=0";
           "x=0.99a, y=0";
           "x=1.01a, y=0";
           "x=2a, y=0";
           "y=0, z=0.01 m";
           "y=0, z=a";
           "y=0, z=2a"];

componentes = ["Bx"; "By"; "Bz"];
fig = 0;

// ---------------- Calculo y 24 graficas ------------------------
for r = 1:8
    if rectas(r,1) == 1 then
        // Rectas paralelas al eje z.
        xp = rectas(r,2)*ones(Np,1);
        yp = zeros(Np,1);
        zp = t;
        variable = t;
        etiqueta_eje = "z (m)";
    else
        // Rectas paralelas al eje x: z permanece fijo.
        xp = t;
        yp = zeros(Np,1);
        zp = rectas(r,2)*ones(Np,1);
        variable = t;
        etiqueta_eje = "x (m)";
    end

    [Bx, By, Bz] = campo_espira(xp, yp, zp, a, I, rectas(r,3), blk);

    // Por simetria respecto del plano xz, By=0 en todas estas rectas.
    // En el eje z (x=y=0), tambien Bx=0 por simetria cilindrica.
    By(:) = 0;
    if rectas(r,1) == 1 & abs(rectas(r,2)) < %eps then
        Bx(:) = 0;
    end

    B = [Bx, By, Bz];

    for c = 1:3
        fig = fig + 1;
        scf(fig);
        plot(variable, B(:,c));
        xgrid();
        xlabel(etiqueta_eje);
        ylabel(componentes(c) + " (T)");
        title(componentes(c) + " - " + nombres(r));
    end
end

mprintf("Se generaron %.0f graficas independientes (%.0f rectas x 3 componentes).\n", fig, 3*8);
mprintf("Resolucion: %.0f puntos por recta, intervalo abierto (-5a,5a).\n", Np);
mprintf("El radio usado es %.3f m y la corriente %.3f A.\n", a, I);
