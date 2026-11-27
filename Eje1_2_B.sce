// ================================================================
// FISICA 2 - Espira circular - Apartado 1.2(b)
// Comparacion en un mismo grafico de:
//   1) Bz analitico sobre el eje de la espira.
//   2) Bz numerico calculado mediante la sumatoria de Biot-Savart.
//
// Sobre el eje: x=0, y=0; intervalo 0 <= z <= 5a.
// Unidades SI: a [m], I [A], z [m], Bz [T].
// ================================================================

clear;
clc;

// ---------------- Parametros fisicos ----------------------------
mu0 = 4*%pi*1e-7;   // permeabilidad del vacio [T*m/A]
a   = 0.03;          // radio [m]
I   = 1.0;           // corriente [A]

// ---------------- Parametros numericos --------------------------
Np = 1000;           // puntos de observacion
M  = 4096;           // divisiones angulares de la integral

// Se incluyen los extremos z=0 y z=5a.
z = linspace(0, 5*a, Np)';

// ---------------- Solucion analitica ----------------------------
// Bz(0,z) = mu0*I*a^2 / (2*(a^2+z^2)^(3/2))
Bz_analitico = (mu0*I*a^2) ./ (2*(a^2 + z.^2).^(3/2));

// ---------------- Solucion numerica -----------------------------
// Ley de Biot-Savart discretizada con regla del punto medio:
// Bz = mu0*I/(4*pi) * sum_j [a^2/(a^2+z^2)^(3/2) * dphi]
// En el eje, cada contribucion angular tiene el mismo valor.
// Se realiza explicitamente la suma de las M contribuciones.
dphi = 2*%pi/M;
factor = mu0*I/(4*%pi)*dphi;
denominador = (a^2 + z.^2).^(3/2);
contribucion = factor * (a^2 ./ denominador);

Bz_numerico = zeros(Np, 1);
for j = 1:M
    Bz_numerico = Bz_numerico + contribucion;
end

// ---------------- Grafica de las dos curvas --------------------
scf(1);
clf();

// Ambas curvas se dibujan en el mismo sistema de ejes:
// Campo teorico: azul continuo.
// Campo numerico (sumatoria): rojo punteado.
// Se usan los 1000 puntos para ambas curvas.
plot(z, Bz_analitico, "b-", z, Bz_numerico, "r:");

xgrid();
xlabel("z (m)");
ylabel("B_z (T)");
title("comparacion entre campos");
legend("Campo teorico", "Campo por sumatoria", "in_upper_right");

// ---------------- Carpeta de exportacion ------------------------
// Busca el Escritorio local o el Escritorio de OneDrive.
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
        if oneDrive == "" then
            oneDrive = getenv("OneDriveConsumer");
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

// Exporta la figura comparativa como PNG.
ruta_png = carpeta_salida + "/Grafica26Bz_comparacion_entre_campos_x=0_y=0.png";
xs2png(1, ruta_png);

// Exporta los 1000 valores para poder revisarlos o procesarlos.
ruta_csv = carpeta_salida + "/Datos_Parte2b_comparacion_Bz.csv";
fid = mopen(ruta_csv, "wt");
if fid == -1 then
    error("No se pudo crear el archivo CSV de resultados.");
end
mfprintf(fid, "z_m;Bz_analitico_T;Bz_numerico_T\n");
for k = 1:Np
    mfprintf(fid, "%.12e;%.12e;%.12e\n", z(k), Bz_analitico(k), Bz_numerico(k));
end
mclose(fid);

// ---------------- Resumen en consola ----------------------------
error_abs_max = max(abs(Bz_numerico - Bz_analitico));

mprintf("\nCOMPARACION ANALITICA Y NUMERICA - APARTADO 1.2(b)\n");
mprintf("Radio a = %.6f m\n", a);
mprintf("Corriente I = %.3f A\n", I);
mprintf("Puntos sobre el eje = %d\n", Np);
mprintf("Divisiones angulares M = %d\n", M);
mprintf("Intervalo: 0 <= z <= %.6f m\n", 5*a);
mprintf("Bz analitico en z=0 = %.12e T\n", Bz_analitico(1));
mprintf("Bz numerico en z=0 = %.12e T\n", Bz_numerico(1));
mprintf("Maxima diferencia absoluta = %.12e T\n", error_abs_max);
mprintf("Grafica guardada en: %s\n", ruta_png);
mprintf("Datos guardados en: %s\n\n", ruta_csv);
