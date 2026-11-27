// ================================================================
// FISICA 2 - Solenoide corto: apartado 2.1(b)
// Comparacion NUMERICA y ANALITICA de Bz(0,0,z) sobre el eje.
//
// Datos del enunciado:
//   Radio: a = 3 cm = 0.03 m
//   Longitud: L = a = 0.03 m
//   Cantidad de espiras: N = 200
//   Corriente: I = 1 A
//   Intervalo: -2L <= z <= 2L
//
// Metodos:
//   Numerico: suma de los campos de las N espiras discretas.
//   Analitico: modelo continuo de solenoide finito (n = N/L).
//
// Convencion: corriente antihoraria vista desde +z.
// Unidades: z en m y Bz en T.
// ================================================================

clear;
clc;

// ---------------- Parametros fisicos ----------------------------
mu0 = 4*%pi*1e-7;  // T*m/A
a   = 0.03;         // m
L   = a;            // m
I   = 1.0;          // A
N   = 200;          // numero de espiras
n   = N/L;          // espiras/m

// ---------------- Parametros numericos --------------------------
Np = 1000;          // cantidad de puntos sobre el eje

dz = L/N;
// Se usa el punto medio de cada tramo de longitud dz.
k = (0:N-1)';
zk = -L/2 + (k+0.5)*dz;

// Intervalo pedido, incluyendo ambos extremos.
z = linspace(-2*L, 2*L, Np)';

// ---------------- Calculo numerico ------------------------------
// Campo axial de una espira centrada en zk:
// Bz_k = mu0*I*a^2 / (2*(a^2+(z-zk)^2)^(3/2))
// Se suman las contribuciones de las 200 espiras.
Bz_num = zeros(Np,1);

for i = 1:Np
    distancia = z(i) - zk;
    contribuciones = mu0*I*a^2 ./ (2*(a^2 + distancia.^2).^(3/2));
    Bz_num(i) = sum(contribuciones);
end

// ---------------- Calculo analitico -----------------------------
// Modelo de un solenoide finito con espiras distribuidas
// de forma continua y uniforme (n = N/L).
Bz_ana = (mu0*n*I/2) * ( ...
    (z+L/2)./sqrt(a^2+(z+L/2).^2) ...
    - (z-L/2)./sqrt(a^2+(z-L/2).^2) ...
);

// ---------------- Comparacion de resultados ---------------------
error_abs = abs(Bz_num - Bz_ana);
error_rel = 100 * error_abs ./ abs(Bz_ana); // porcentaje

// ---------------- Grafica ---------------------------------------
scf(1);
clf();
plot(z, [Bz_num Bz_ana]);
xgrid();
xlabel("z (m)");
ylabel("Bz (T)");
title("Apartado 2.1(b) - Solenoide corto: numerico vs analitico");
legend(["Numerico (200 espiras)", "Analitico (solenoide finito)"], "in_upper_right");

// ---------------- Carpeta de exportacion ------------------------
// Busca el Escritorio de Windows (incluyendo OneDrive).
// Como alternativa, prueba ~/Desktop o ~/Escritorio en Linux/macOS.
perfil_usuario = getenv("USERPROFILE");
escritorio = "";

if perfil_usuario <> "" then
    if isdir(perfil_usuario + "/Desktop") then
        escritorio = perfil_usuario + "/Desktop";
    elseif isdir(perfil_usuario + "/Escritorio") then
        escritorio = perfil_usuario + "/Escritorio";
    end

    if escritorio == "" then
        oneDrive = getenv("OneDrive");
        if oneDrive == "" then
            oneDrive = getenv("OneDriveCommercial");
        end
        if oneDrive <> "" then
            if isdir(oneDrive + "/Desktop") then
                escritorio = oneDrive + "/Desktop";
            elseif isdir(oneDrive + "/Escritorio") then
                escritorio = oneDrive + "/Escritorio";
            end
        end
    end
else
    home_usuario = getenv("HOME");
    if home_usuario <> "" then
        if isdir(home_usuario + "/Desktop") then
            escritorio = home_usuario + "/Desktop";
        elseif isdir(home_usuario + "/Escritorio") then
            escritorio = home_usuario + "/Escritorio";
        end
    end
end

if escritorio == "" then
    error("No se encontro el Escritorio. Edita la variable escritorio e indica su ruta.");
end

carpeta_salida = escritorio + "/Graficas_Fisica";
if ~isdir(carpeta_salida) then
    [estado_mkdir, mensaje_mkdir] = mkdir(carpeta_salida);
    if estado_mkdir <> 1 & estado_mkdir <> 2 then
        error("No se pudo crear la carpeta de salida: " + mensaje_mkdir);
    end
end

// Guarda la grafica PNG.
ruta_png = carpeta_salida + "/Solenoide_2_1b_comparacion.png";
xs2png(1, ruta_png);

// Guarda los resultados en CSV para Excel.
ruta_datos = carpeta_salida + "/Datos_Solenoide_2_1b.csv";
fid = mopen(ruta_datos, "wt");
if fid == -1 then
    error("No se pudo crear el archivo de datos: " + ruta_datos);
end
mfprintf(fid, "z_m;Bz_numerico_T;Bz_analitico_T;error_relativo_porcentaje\n");
for j = 1:Np
    mfprintf(fid, "%.12e;%.12e;%.12e;%.12e\n", ...
        z(j), Bz_num(j), Bz_ana(j), error_rel(j));
end
mclose(fid);

// ---------------- Resumen en consola ----------------------------
// Como Np es par, z=0 no cae exactamente en la malla;
// calculamos Bz en z=0 por separado.
Bz_centro_num = sum(mu0*I*a^2 ./ (2*(a^2+zk.^2).^(3/2)));
Bz_centro_ana = (mu0*n*I/2) * ( ...
    (L/2)/sqrt(a^2+(L/2)^2) - (-L/2)/sqrt(a^2+(-L/2)^2) ...
);

mprintf("\nAPARTADO 2.1(b) - COMPARACION NUMERICA Y ANALITICA\n");
mprintf("Radio a = %.6f m\n", a);
mprintf("Longitud L = %.6f m\n", L);
mprintf("Corriente I = %.3f A\n", I);
mprintf("Espiras N = %d\n", N);
mprintf("Densidad n = %.6f espiras/m\n", n);
mprintf("Puntos sobre el eje = %d\n", Np);
mprintf("Intervalo: %.6f <= z <= %.6f m\n", -2*L, 2*L);
mprintf("Bz numerico en z=0 = %.12e T\n", Bz_centro_num);
mprintf("Bz analitico en z=0 = %.12e T\n", Bz_centro_ana);
mprintf("Error relativo maximo en la malla = %.8e %%\n", max(error_rel));
mprintf("Grafica guardada en: %s\n", ruta_png);
mprintf("Datos guardados en: %s\n\n", ruta_datos);
