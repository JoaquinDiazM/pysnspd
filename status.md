# Estado de pySNSPD 1.1.1

**Versión vigente 1.1.1:** rigidez dependiente de Te, D calibrada y calor firmado compatible con Pauli. Dos pruebas nuevas de 0.2 ps completas; bloqueo eliminado, sin admisión de latencia completa. [Diagnóstico y límites](docs/releases/v1.1.1.md).

## Evidencia histórica de v1.0.0

**Etapa 5 completada: admisión exploratoria del modelo efectivo, 12 de octubre de 2026.** La entrada pública es `pysnspd.effective`; el modelo de la memoria permanece congelado en thesis-`v1.0.0`.

| Etapa | Estado final |
|---|---|
| 5.1 | Un único sistema acoplado: gap complejo, ocupaciones electrónicas/fonónicas, potencial y circuito CM |
| 5.2 | Energía común, potencia independiente, Joule/condensado/escape/cortes y ledger integrado auditados |
| 5.3 | Controles temporal, malla, dominio y espectro de 775 nm admitidos en 0–36,6 ps |

La decisión científica explícita está en `/home/jdiaz/scratch/stage5_final_evidence_completed_20261012/admission.json`, SHA-256 `ed30c8525ecb3952e4042eca3b84e9ea0c2c78f41d7e0f17d880779b68747c4a`. La revisión manual de cobertura tiene SHA-256 `2c404823c78594ae22fad71f1fee50e62b5d9ad3ee7beab89436b9bfd40252ac`. Estos certificados identifican originales, fuentes, parámetros, decisiones y evidencia; no se infiere admisión sólo de un PASS automático.

## Resultados nominales

| Caso | Ventana y resultado | Máximo defecto energético global |
|---|---|---|
| 775 nm | Hotbelt guardada a 15,3 ps; cruce 31,94131686165 ps; confirmación cruce+0,5 ps; final 36,94131686165 ps | 0,663463 % |
| 1550 nm | 200 ps sin cruce de 0,1 mV; latencia sin asignar; sin hotbelt transversal según el criterio guardado | 0,000144050 % |

El defecto se normaliza por $E_{\rm ret}+\sum|\Delta W_{\rm ext}|$ y debe ser ≤1 %. Las auditorías individuales pasan energía instantánea/integrada, fuente única, soporte Pauli/Bose, estados/eventos guardados y corrientes globales positivas. La señal basal DC satisface el criterio nV. La hotbelt es la componente conectada con $|\Delta|\le0.2|\Delta_{\rm DC}|$ que toca ambos bordes transversales; es un diagnóstico separado del gatillo eléctrico.

## Comparaciones numéricas finales: 775 nm, 0–36,6 ps

| Eje | Error Vout [%] | L2 de campo/DC [%] | Déficit/pico real [%] | Diferencia de latencia [ps] | Resultado |
|---|---:|---:|---:|---:|---|
| Temporal | 0.089801 | 0.149759 | 0.053830 | 0.007566 | PASS |
| Malla | 4.039153 | 4.817296 | 3.015699 | 0.517963 | PASS |
| Dominio | 0.843195 | 2.737720 | 1.256287 | 0.115053 | PASS |
| Espectro | 2.597424 | 3.085255 | 1.971201 | 0.313528 | PASS |

El temporal exige Vout y mínimo nodal/DC al 2 % (este último: 0,639400 %); su latencia usa máximo entre 2 % y 0,1 ps. Malla y dominio usan 5 % en Vout, latencia y ambos observables espaciales. En espectro, campo/DC y déficit/pico usan 5 %; Vout usa el 5 % adoptado explícitamente **después de los datos**, mientras latencia conserva 2 %. El mínimo nodal es diagnóstico en malla/dominio/espectro; el L2 temporal mostrado es diagnóstico adicional.

El FAIL espectral previo de Vout 2,597424 % frente al 2 % permanece íntegro. También permanecen los resultados originales 2 % y 3/2/2, los máximos históricos y el diagnóstico nodal de malla 22,5717 %/DC. No se declara precisión asintótica por adoptar estos márgenes de ingeniería.

La revisión cubre formación y confirmaciones hasta 33 ps, con al menos 3,640719 ps comparados tras la última confirmación. **Todas las colas posteriores a 36,6 ps quedan fuera de la convergencia de cuatro ejes**, aunque cada original completa cruce+5 ps y conserva auditoría individual. No se certifica el pico completo ni recuperación/reset total.

## Modelo y supuestos admitidos

Referencia: ancho 80 nm, espesor físico 7 nm, separación de malla 5 nm, Tb=0,9 K e I=21,297206 µA. `h_mesh` es resolución, no espesor. La producción exploratoria utiliza **54 estados electrónicos**; 80 estados son un control opcional de sensibilidad, no verdad continua.

El cierre neutral BCS/KWT efectivo conserva su energía común, usa rigidez radial fenomenológica y material fonónico acústico. τ_ee(Tc)=5 ps y τ_ep(Tc)=24,7 ps son nominales libres de Allmaras (2020), no valores fijos ni centros de un rango requerido. η=0,667, σ=7,549834435 nm y τ_escape=15 ps son hipótesis exploratorias. El estado inicial puede tener temperatura equivalente superior a Tc: KWT y rigidez evaluada al baño son extrapolaciones fenomenológicas calibrables.

El reloj empieza en la transferencia **postcascada**; no se añade 0,57 ps a la latencia. La fuente FD/BE común distribuye la energía entre electrones y fonones. Los controles σ/√2 y √2σ completaron 20 ps y cambian tanto extensión como reparto inicial: no son una prueba geométrica pura ni censuras a 200 ps.

1550 nm tiene una trayectoria nominal conservativa y contraste temporal tardío propio 29,62–31,62 ps; no posee cuatro gemelos refinados de 0–200 ps. No se afirma jitter estadístico, latencia óptica absoluta, reset eléctrico total, reproducción experimental calibrada ni tasas absolutas deducidas del archivo MIT. El intento Moyal rechazado y los DC/controles previos conservan sus decisiones históricas.

## Reproducción e historia

[Modelo vigente](docs/implementation/MODELO_VIGENTE.md), [etapa 5](docs/implementation/stage5/README.md), [comandos normales/visuales](docs/GEMINGA_COMMANDS_PRESENT.md) e [informe final](output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf). Las ecuaciones están en [Principal](Matematica/Principal.pdf), [Auxiliar](Matematica/Auxiliar.pdf) y [Precomputo](Matematica/Precomputo.pdf).

La [versión anterior de este estado](docs/releases/history/status_before_v1.1.0_20261011.md) permanece histórica. El cuaderno GEMINGA_COMMANDS.md conserva comandos posibles y usados. Datos pesados, builds, cachés y renders permanecen en `/home/jdiaz/scratch`; se preservan ambos enlaces tmp. El presupuesto global de ejecución no supera el 90 % disponible.
