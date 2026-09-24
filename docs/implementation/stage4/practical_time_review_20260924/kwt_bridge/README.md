# Reutilización directa del paso KWT de la memoria

**Comprobación completada; producción intacta.** El adaptador experimental
`pysnspd/experimental/heredado_kwt_bridge.py` llama directamente a
`TDGLSolver.solve_for_psi_squared`. Conserva las unidades y la movilidad de
`ThermalKWTNormal`, mediante P.1–P.2 de
[la revisión de métodos publicados](../published_methods.md).

La entrada `gradient` es la fuerza integrada por nodo de la nueva acción.
`current` es su corriente espectral orientada por arista. El adaptador obtiene
el potencial con el bloque normal del mismo modelo y convierte el paso físico
en picosegundos al tiempo del backend heredado. El argumento obligatorio
`delta0_over_kBTc` indica la normalización fija del gap; no modifica la amplitud
física. No se añade corriente GL ni se aplica dos veces la movilidad.

El control usa un campo suave, con amplitud entre 1,55 y 1,63 en unidades
`kB Tc` y una fase suave de hasta 0,12 radianes, sobre 63 nodos. Los 28 nodos
del contorno están fijos. Las 16 raíces espectrales se obtienen con tolerancia
10⁻⁷; su máximo residuo fue 9,52×10⁻⁸. Esta suma es una entrada para verificar
la conversión temporal, no un nuevo cálculo convergido de propiedades NbN.

Se compara `(d_nuevo-d_actual)/dt_ps` con la velocidad independiente
`ThermalKWTNormal.response(...)[velocity]/tD_ps`. La norma es
`sqrt(sum_i m_i |error_i|²)`, con pesos de área adimensionales `m_i`; el error
relativo divide por la norma de la velocidad completa de referencia.

| Paso físico (ps) | Error relativo de la velocidad | Cambio de los contactos |
|--:|--:|--:|
| 0,0001 | 0,01878 % | 0 exacto |
| 0,00005 | 0,009389 % | 0 exacto |
| 0,000025 | 0,004697 % | 0 exacto |

Al dividir el paso por dos, el error se reduce por factores **2,00003** y
**1,99893**, coherentes con el límite de primer orden del Euler heredado. La
fuerza y el potencial son distintos de cero: no se usa un equilibrio que
ocultaría un factor incorrecto. El potencial máximo es 0,3963 en unidades
`v = 2 e phi/(kB Tc)`. La contribución del enlace temporal a la norma de la
velocidad es 1,758 ps⁻¹, frente a 17,81 ps⁻¹ de la movilidad material; ambas
intervienen efectivamente en el contraste.

El [recibo reproducible](receipt.json) conserva cifras, versiones, tiempos y
SHA-256 de los módulos utilizados. El cálculo demoró 0,172 s. El test focalizado
`tests/test_heredado_kwt_bridge.py` pasó en 0,181 s y verifica las tres llamadas
al método heredado real.

Ambas comprobaciones se ejecutaron localmente en Windows, con Python 3.12.14,
NumPy 2.5.3 y SciPy 1.18.1, limitando las bibliotecas a un hilo. El recibo conserva
los hashes exactos de ese checkout. `pysnspd/solver/core.py` tiene allí saltos de
línea CRLF; su contenido normalizado coincide byte a byte con Git y con el
archivo LF de Geminga. No se alteró el método heredado.

Para repetir el diagnóstico en Geminga con un archivo nuevo:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -m sandbox.stage4_core.kwt_bridge_check \
  --output tmp/kwt_bridge_repeat_20260924.json
```

La prueba admite reutilizar **el paso local** con esta conversión. No compara
una trayectoria completa Euler–ETD2, no decide el paso máximo del dispositivo
y no integra todavía transporte, puertos o circuito. La corriente futura debe
continuar viniendo de la nueva acción y conectarse con sus factores SI. Cambiar
solo la fuerza y dejar activa la corriente GL del solver completo no está
admitido por este puente.
