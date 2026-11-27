// ================================================================
// FISICA 2 - Espira circular - Apartado 1.2(c)
// Error relativo entre Bz numerico (sumatoria de Biot-Savart)
// y Bz analitico sobre el eje de la espira.
//
// Punto de observacion: x=0, y=0.
// Intervalo: 0 <= z <= 5a.
// Se grafica el error relativo en funcion de z para varias
// resoluciones angulares y se cuantifica el error maximo y medio.
// Unidades SI: a [m], I [A], z [m], Bz [T].
// ================================================================

clear;
clc;

// ---------------- Parametros fisicos ----------------------------
mu0 = 4*%pi*1e-7;    // permeabilidad del vacio [T*m/A]
a   = 0.03;           // radio de la espira [m]
I   = 1.0;            // corriente [A]

// ---------------- Parametros numericos --------------------------
Np = 1000;                            // puntos sobre el eje
M_values = [256 512 1024 2048 4096 8192]; // resoluciones angulares
M_referencia = 4096;                  // resolucion usada en la grafica principal

// Incluye los extremos z=0 y z=5a.
z = linspace(0, 5*a, Np)';

// ---------------- Campo analitico -------------------------------
// Bz(0,z) = mu0*I*a^2 / (2*(a^2+z^2)^(3/2))
denominador = (a^2 + z.^2).^(3/2);
Bz_analitico = (mu0*I*a^2) ./ (2*denominador);

// ---------------- Sumatoria numerica y errores ------------------
nM = length(M_values);
errores_pct = zeros(Np, nM);
error_max_pct = zeros(1, nM);
error_medio_pct = zeros(1, nM);
Bz_numerico_referencia = zeros(Np, 1);

for r = 1:nM
    M = M_values(r);
    dphi = 2*%pi/M;

    // En el eje, x=y=0, la componente z del integrando es:
    // a^2/(a^2+z^2)^(3/2), independiente del angulo phi.
    // Se conserva la sumatoria explicita de M contribuciones.
    factor = mu0*I/(4*%pi)*dphi;
    contribucion = factor*(a^2 ./ denominador);
    Bz_numerico = zeros(Np, 1);

    for j = 1:M
        Bz_numerico = Bz_numerico + contribucion;
    end

    // Error relativo expresado en porcentaje.
    // error_relativo(%) = 100*|Bnum-Banal|/|Banal|
    error_relativo_pct = 100*abs(Bz_numerico - Bz_analitico) ./ abs(Bz_analitico);

    errores_pct(:, r) = error_relativo_pct;
    error_max_pct(r) = max(error_relativo_pct);
    error_medio_pct(r) = mean(error_relativo_pct);

    if M == M_referencia then
        Bz_numerico_referencia = Bz_numerico;
    end
end

// ---------------- Grafica 1: error relativo vs z ----------------
// Sobre el eje, la sumatoria angular es exactamente igual a la
// expresion analitica en aritmetica ideal. La computadora puede
// producir diferencias diminutas por redondeo (del orden de 1e-11 %).
// Para que la escala de la grafica no amplifique ese ruido, los valores
// inferiores a 1e-10 % se muestran como cero SOLO en la grafica.
// Los valores numericos originales, sin modificar, se conservan en CSV
// y se informan en la consola.
umbral_grafica_pct = 1e-10;
indice_referencia = find(M_values == M_referencia);
error_grafica_pct = errores_pct(:, indice_referencia);

for k = 1:Np
    if abs(error_grafica_pct(k)) < umbral_grafica_pct then
        error_grafica_pct(k) = 0;
    end
end

scf(1);
clf();
plot(z, error_grafica_pct, "b-");
ax = gca();
ax.data_bounds = [0, -umbral_grafica_pct; 5*a, umbral_grafica_pct];
xgrid();
xlabel("z (m)");
ylabel("Error relativo (%)");
title("Error relativo sobre el eje: practicamente E_r = 0 (M = 4096)");
legend("Error relativo numerico (valores < 1e-10 % mostrados como 0)", "in_upper_right");

// ---------------- Grafica 2: error maximo vs M ------------------
// Muestra el pequeno error de redondeo REAL calculado, sin umbralizar.
// No es error de integracion: el integrando angular es constante.
scf(2);
clf();
plot(M_values, error_max_pct, "bo-");
ax = gca();
ax.auto_clear = "off";
plot(M_values, error_medio_pct, "rs-");
ax = gca();
ax.auto_clear = "on";
xgrid();
xlabel("Resolucion angular M (divisiones)");
ylabel("Error relativo (%)");
title("Error de redondeo maximo y medio segun la resolucion angular");
legend("Error maximo real", "Error medio real", "in_upper_right");

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
    error("No se encontro el Escritorio. Edita la variable escritorio e indica la ruta.");
end

carpeta_salida = escritorio + "/Graficas_Fisica";
if ~isdir(carpeta_salida) then
    [estado_mkdir, mensaje_mkdir] = mkdir(carpeta_salida);
    if estado_mkdir <> 1 & estado_mkdir <> 2 then
        error("No se pudo crear la carpeta de salida: " + mensaje_mkdir);
    end
end

// Exporta las dos graficas.
ruta_png_error = carpeta_salida + "/Grafica27Bz_error_relativo_vs_z.png";
xs2png(1, ruta_png_error);

ruta_png_convergencia = carpeta_salida + "/Grafica28Bz_error_vs_resolucion.png";
xs2png(2, ruta_png_convergencia);

// ---------------- Exportacion de resultados CSV -----------------
// Archivo con los errores relativos en cada z para todas las M.
ruta_csv = carpeta_salida + "/Datos_Parte2c_error_relativo.csv";
fid = mopen(ruta_csv, "wt");
if fid == -1 then
    error("No se pudo crear el archivo CSV de errores relativos.");
end

mfprintf(fid, "z_m;Bz_analitico_T;Bz_numerico_M4096_T");
for r = 1:nM
    mfprintf(fid, ";error_relativo_M%d_pct", M_values(r));
end
mfprintf(fid, "\n");

for k = 1:Np
    mfprintf(fid, "%.12e;%.12e;%.12e", z(k), Bz_analitico(k), Bz_numerico_referencia(k));
    for r = 1:nM
        mfprintf(fid, ";%.12e", errores_pct(k,r));
    end
    mfprintf(fid, "\n");
end
mclose(fid);

// Archivo resumen de la prueba de resolucion.
ruta_csv_conv = carpeta_salida + "/Datos_Parte2c_convergencia.csv";
fid = mopen(ruta_csv_conv, "wt");
if fid == -1 then
    error("No se pudo crear el archivo CSV de convergencia.");
end
mfprintf(fid, "M;error_maximo_pct;error_medio_pct\n");
for r = 1:nM
    mfprintf(fid, "%d;%.12e;%.12e\n", M_values(r), error_max_pct(r), error_medio_pct(r));
end
mclose(fid);

// ---------------- Resumen en consola ---------------------------
mprintf("\nPARTE 1.2(c) - ERROR RELATIVO Y RESOLUCION NUMERICA\n");
mprintf("En la grafica 27, errores menores que %.1e %% se muestran como cero para visualizar E_r ~= 0.\n", umbral_grafica_pct);
mprintf("Los CSV y esta tabla conservan los errores reales sin modificar.\n");
mprintf("Radio a = %.6f m\n", a);
mprintf("Corriente I = %.3f A\n", I);
mprintf("Puntos sobre el eje = %d\n", Np);
mprintf("Intervalo: 0 <= z <= %.6f m\n", 5*a);
mprintf("\nResolucion M | Error maximo (%%) | Error medio (%%)\n");
for r = 1:nM
    mprintf("%10d | %16.8e | %15.8e\n", M_values(r), error_max_pct(r), error_medio_pct(r));
end
mprintf("\nResolucion de referencia para la curva principal: M = %d\n", M_referencia);
mprintf("Error real maximo para M=%d: %.12e %%\n", M_referencia, error_max_pct(indice_referencia));
mprintf("Grafica de error guardada en: %s\n", ruta_png_error);
mprintf("Grafica de resolucion guardada en: %s\n", ruta_png_convergencia);
mprintf("Datos detallados guardados en: %s\n", ruta_csv);
mprintf("Resumen de convergencia guardado en: %s\n\n", ruta_csv_conv);
