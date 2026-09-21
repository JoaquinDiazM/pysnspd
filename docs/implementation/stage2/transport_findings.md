# Transporte entre dos celdas

El transporte satisface el presupuesto de 0,1 % con la malla complementaria de
630 estados por celda en los seis perfiles registrados. El catálogo R2 permanece
intacto. La malla complementaria obtiene la energía y sus derivadas del mismo
conteo causal de estados; no interpola la fuerza de manera independiente.

| Estados por celda | Error máximo de respuesta | Error máximo de potencia |
|---:|---:|---:|
| 630 | 0,0039463 % | 0,00200317 % |
| 1260 | 0,000990867 % | 0,000501739 % |
| 2520 | 0,000248253 % | 0,000125556 % |

Los errores disminuyen aproximadamente cuatro veces al duplicar la resolución.
Las referencias independientes de cuadratura 24 y 48 coinciden a menos de
7,05·10⁻¹⁶ relativo. La prueba completa tardó 170,26 s en Geminga; no integra una
trayectoria temporal y ese tiempo no se extrapola al costo de un transitorio.

La conservación de número de cuasipartículas y energía, el equilibrio de Fermi
común a las dos celdas, las fronteras de Pauli y la producción de entropía pasan
sus controles en las tres resoluciones. La mayor respuesta residual ponderada
en equilibrio es 6,85·10⁻¹⁸.

Se conserva el resultado inicial fallido: con 180 estados, el perfil no térmico
con brechas distintas daba 0,25534 % de error de respuesta y 0,137989 % de potencia.
Subir el orden de integración de 2 a 4 apenas lo cambiaba. Una descomposición
adicional mostró que retirar el factor de atenuación geométrico tampoco resolvía
el error: dominaba la representación del perfil entre estados. No se modificó
la ley de transporte ni se normalizaron sus productos de Pauli.

Los casos incluyen estado normal, superconductividad con brecha y sin brecha;
sus parámetros y perfiles térmicos/no térmicos están registrados en
[transport_results_scenarios.json](transport_results_scenarios.json). Los números
y hashes están en [transport_results.json](transport_results.json), y el fallo
inicial en [pilot_initial/transport_original_results.json](pilot_initial/transport_original_results.json).
Los ensayos complementarios anteriores permanecen en
[history/transport_edge2](history/transport_edge2/snapshot.json) y
[history/transport_edge8](history/transport_edge8/snapshot.json). La cuadratura
final también refina el borde y el intervalo de conteo 0,01–0,1, exigidos por
controles independientes de fuerzas y perfiles estrechos.

La resolución para la dinámica simultánea debe satisfacer también las pruebas
de colisiones. Estos ensayos no certifican cualquier distribución arbitrariamente
estrecha, transporte fuera del soporte energético ni un transitorio completo.
