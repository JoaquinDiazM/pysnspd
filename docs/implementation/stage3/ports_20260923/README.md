# Etapa 3: cierre estático nodal y entrada de bordes/circuito

Resultado vigente del 23 de septiembre de 2026. [Informe con plots](Informe_avance_bordes_etapa_3_20260923.md)
([PDF](../../../../output/pdf/implementation/Informe_avance_bordes_etapa_3_20260923.pdf)).

- 3A: cerrado para desarrollo en los estados estáticos registrados. Campaña 18/18,
  diez estimaciones relativas aceptadas, una sin certificado relativo, ningún fallo.
  Mayor estimación espacial 0,613214% frente al objetivo 1%.
- Nuevos módulos: funcional abierto 1D, rama térmica de reservorio, intercambio
  electrónico a energía compartida, potencial conservativo y circuito de memoria.
- Verificación nueva en Geminga: 168 pruebas y 23 subpruebas; 9,87 s.
- Piloto abierto: 17 nodos/360 nm, 20,9 s, PASS; error corriente 0,0006035%,
  gradiente terminal 0,07130%. La hélice inicial ya cumple: no demuestra relajación.
- Control circuital independiente: resistencia prescrita 1000 ohm, inductancia
  exterior explícita 7 nH; error Vout 6,45e-14 V frente a exponencial matricial.
- Pendiente manual: seis casos abiertos, tres longitudes y dos resoluciones,
  registrados antes de ver los resultados; estimación 6-10 min. Comando vigente en
  `/home/jdiaz/GEMINGA_COMMANDS.md`. No fue lanzado por el agente.

La etapa 3 completa, interfaz 2D-1D, transiente débil acoplado y producción siguen
pendientes. No se atribuye una señal SNSPD al circuito con resistencia prescrita.
La partición 10 nH = 0,166621 nH resueltos + 9,833379 nH exteriores corresponde al
piloto uniforme 360 nm a 8,63351 microA. No es la inductancia 7 nH del control aislado.

## Evidencia y reproducción

`nodal_audit.json` y `nodal_closure.md` auditan la campaña completa preservada en
`raw/nodal_campaign`. `open_registration.json` fija el siguiente ensayo;
`open_pilot` conserva fuentes, resultados, rama, iteraciones y progreso.
`electrical` incluye datos, figuras y balances del control de puerto/circuito.
`tests.log` registra la suite nueva. `validation_review.json` es una revisión
independiente de archivos y aritmética guardada, sin cálculo físico nuevo.

El funcional abierto usa elementos de Lobatto de grado 4, cuadratura positiva y
operadores compatibles. Mantiene un espectro por nodo y una energía común para
fuerza/corriente. La energía libre térmica incluye entropía; no se diferencia la
energía interna recalculando FD como si fuera una variación a ocupación fija.
Los extremos fijan amplitud de reservorio e imponen corriente por trabajo
variacional. Fijar fase derecha 0 sólo elimina el gauge global, no impone una
diferencia de fase física adicional. El potencial tiene su propio gauge eléctrico.

El intercambio electrónico con el baño fija una FD a la misma temperatura y
evalúa eventos a energía física común. Su geometría y escala temporal son
explícitas; el baño tiene un balance opuesto al de la celda, no un volumen
finito inventado. Las pruebas no garantizan positividad de un integrador arbitrario.

## Alcance del próximo lote

Compara 360/720/1080 nm con 4/8, 8/16, 12/24 elementos de grado 4 respectivamente.
Las prolongaciones son 0/180/360 nm por extremo del segmento central. Evalúa
equilibrio, flujo terminal, soporte y signo D.36; no prueba reflexión dinámica,
acoplamiento 2D-1D ni evolución del detector. Conserva salidas sin sobrescribir
y se detiene ante fallo sin cambiar solver, tolerancias o estado.

Las libretas y entradas previas se conservaron exactas en
`GEMINGA_COMMANDS_before_ports.md` y `previous_delivery_exact/`. El verificador
actual incluye el manifiesto previo, con esas dos rutas históricas resueltas
a su instantánea. No se alteró la entrega inmutable de etapa 2 ni el tag v1.0.0.
