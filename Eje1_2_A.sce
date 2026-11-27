// ================================================================
// FISICA 2 - Espira circular: apartado 1.2(a)
// Calculo NUMERICO de Bz(0,0,z) sobre el eje de la espira.
//
// Datos del enunciado:
//   Radio: a = 3 cm = 0.03 m
//   Corriente: I = 1 A
//   Intervalo: 0 <= z <= 5a
//
// Metodo:
//   Se aplica la ley de Biot-Savart y se discretiza la integral
//   sobre el angulo phi mediante la regla del punto medio.
//   La integral se aproxima con M contribuciones angulares.
//
// Convencion: corriente antihoraria vista desde +z; en el eje,
// el campo apunta en la direccion +z (Bz positivo).
// Unidades: z en m y Bz en T.
// ================================================================

clear;
clc;

// ---------------- Parametros fisicos ----------------------------
mu0 = 4*%pi*1e-7;  // T*m/A
a   = 0.03;         // m
I   = 1.0;          // A

// ---------------- Parametros numericos --------------------------
Np  = 1000;         // cantidad de puntos de observacion sobre el eje
M   = 4096;         // divisiones angulares para la integral de Biot-Savart
blk = 50;           // puntos de observacion procesados por bloque

// El enunciado pide incluir ambos extremos: 0 y 5a.
z = linspace(0, 5*a, Np)';

// ---------------- Funcion de calculo numerico -------------------
// En el eje x=0, y=0, la componente Bz de cada elemento es
//   a^2 / (a^2 + z^2)^(3/2) dphi.
// Se aproxima la integral mediante M contribuciones angulares iguales.
function Bz = campo_Bz_eje(zp, a, I, M, blk)
    mu0 = 4*%pi*1e-7;
    n = length(zp);
    dphi = 2*%pi/M;

    // En el eje x=0, y=0, el integrando de Bz es independiente
    // del angulo phi: cada una de las M contribuciones es identica.
    // Por eso, la suma discreta de Biot-Savart se puede evaluar como
    // M veces una contribucion, sin crear matrices angulares.
    Bz = zeros(n,1);
    factor = mu0*I/(4*%pi)*dphi;

    for i = 1:blk:n
        j = min(i+blk-1, n);
        zi = matrix(zp(i:j), j-i+1, 1);

        // Para x=y=0: |r-r'|^2 = a^2 + z^2 y el numerador de Bz es a^2.
        denominador = (a^2 + zi.^2).^1.5;
        Bz(i:j) = factor*M*a^2 ./ denominador;
    end
endfunction

// ---------------- Calculo numerico ------------------------------
Bz = campo_Bz_eje(z, a, I, M, blk);

// ---------------- Grafica ---------------------------------------
scf(1);
plot(z, Bz);
xgrid();
xlabel("z (m)");
ylabel("Bz (T)");
title("Apartado 1.2(a) - Campo numerico Bz sobre el eje de la espira");

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
ruta_png = carpeta_salida + "/Grafica25Bz_x=0_y=0_eje_0a5a.png";
xs2png(1, ruta_png);

// Guarda tambien los 1000 valores numericos en un archivo de texto
// separado por punto y coma, que se puede abrir con Excel.
ruta_datos = carpeta_salida + "/Datos_Parte2a_Bz_eje.csv";
fid = mopen(ruta_datos, "wt");
if fid == -1 then
    error("No se pudo crear el archivo de datos: " + ruta_datos);
end
mfprintf(fid, "z_m;Bz_T\n");
for k = 1:Np
    mfprintf(fid, "%.12e;%.12e\n", z(k), Bz(k));
end
mclose(fid);

// ---------------- Resumen en consola ----------------------------
Bz_centro = Bz(1);
mprintf("\nAPARTADO 1.2(a) - CALCULO NUMERICO\n");
mprintf("Radio a = %.6f m\n", a);
mprintf("Corriente I = %.3f A\n", I);
mprintf("Puntos sobre el eje = %d\n", Np);
mprintf("Divisiones angulares M = %d\n", M);
mprintf("Intervalo: 0 <= z <= %.6f m\n", 5*a);
mprintf("Bz numerico en z=0 = %.12e T\n", Bz_centro);
mprintf("Grafica guardada en: %s\n", ruta_png);
mprintf("Datos guardados en: %s\n\n", ruta_datos);
