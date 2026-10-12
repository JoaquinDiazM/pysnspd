# Etapa 5 completada: modelo efectivo neutral BCS/KWT

**Admisión exploratoria del 12 de octubre de 2026.** Se cierran 5.1 (evolución conjunta), 5.2 (balance independiente) y 5.3 (controles numéricos de 775 nm). La entrada pública es [`pysnspd.effective`](../../../pysnspd/effective.py); [MODELO_VIGENTE](../MODELO_VIGENTE.md) describe ecuaciones, supuestos y parámetros.

## Decisión y trazabilidad

Certificado autoritativo: `/home/jdiaz/scratch/stage5_final_evidence_completed_20261012/admission.json`, SHA-256 `ed30c8525ecb3952e4042eca3b84e9ea0c2c78f41d7e0f17d880779b68747c4a`.

Cobertura manual: `/home/jdiaz/scratch/stage5_final_evidence_completed_20261012/finite_control_coverage_review.json`, SHA-256 `2c404823c78594ae22fad71f1fee50e62b5d9ad3ee7beab89436b9bfd40252ac`. El certificado identifica el auditor independiente de primarios, las cuatro ventanas, las fuentes congeladas y el protocolo real de guardas/recursos. Los 13 pins físicos v6 permanecen congelados; los arreglos posteriores afectan análisis/procedencia/visualización, no el RHS.

La resolución por defecto es **54 estados electrónicos exploratorios** (corte/orden ξ 32/6, fonones 96, transporte 48). El contraste con 80 estados (64/8, 144, 96) es un control opcional de sensibilidad más costoso; no es verdad continua ni traslada automáticamente los cuatro controles a otros parámetros.

## Qué se comprobó

- Un único sistema acopla Δ, p, n, neutralidad eléctrica y los tres estados del circuito CM.
- El calor de condensado y Joule se obtiene físicamente, sin usar el residuo como calor. Colisiones/transporte comparten estequiometría/energía; se registran escape y cortes.
- Primarios y controles completos pasan soporte, corriente positiva, energía instantánea/integrada y eventos guardados. El margen global energético es 1 %.
- Los cuatro ejes 775 nm pasan en **0–36,6 ps**, con cobertura crítica manual hasta **33 ps**.

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

## Resultados físicos y reloj

775 nm presenta hotbelt transversal por primera vez guardada a 15,3 ps (bracket 15,2–15,3), cruza Vout=0,1 mV a 31,94131686165 ps y termina a 36,94131686165 ps. La confirmación ocurre a cruce+0,5 ps; los cinco ps adicionales se cuentan **desde el cruce**. La hotbelt es una componente conectada con gap ≤0,2 de su DC local que toca ambos bordes; es un diagnóstico diferente del evento eléctrico.

1550 nm completa 200 ps sin cruce: su latencia queda sin asignar. Es una trayectoria nominal exploratoria conservativa, con contraste temporal tardío propio 29,62–31,62 ps; no tiene certificación de cuatro ejes durante 0–200 ps. No se le atribuye por extrapolación una eficiencia experimental.

Los defectos globales máximos de primarios son 0,663463 % (775) y 0,000144050 % (1550), respecto a energía retenida más trabajo externo absoluto acumulado. El reloj comienza en la transferencia efectiva postcascada; no se suma 0,57 ps como retraso óptico calculado.

Los ensayos de anchura σ/√2 y √2σ completaron 20 ps con igual energía retenida. FD/BE a temperatura común cambia también el reparto electrónico/fonónico inicial, además de la extensión. Por ello el pico inicial no ordena por sí solo Vout/hotbelt. Son sensibilidades de preparación efectiva en 0–20 ps; no censuras a 200 ps ni mediciones de latencia.

## Supuestos y utilización

Ancho 80 nm; **espesor físico 7 nm**; `h_mesh=5/3,5 nm` sólo modifica resolución. τ_ee(Tc)=5 ps y τ_ep(Tc)=24,7 ps son nominales calibrables de Allmaras (2020), no valores fijos ni centros impuestos. Material acústico, rigidez radial/KWT fuera del equilibrio, σ, η y τ_escape son hipótesis explícitas para exploración y calibración posterior. La fuente MIT externa se conserva documentalmente; no determina las tasas absolutas de esta versión. No hubo cálculo DFPT propio.

El cierre es BCS/KWT **efectivo neutral**, no una evolución microscópica espacial Usadel/Moyal exacta. La temperatura equivalente inicial puede superar Tc; la rigidez evaluada al baño y la movilidad se usan como cierre fenomenológico. Se preservan el circuito completo, sus constantes de tiempo y el signo global de corriente; no hay reset total en los comandos presentes.

El modelo queda disponible para barridos de diseño y comparación experimental bajo hipótesis explícitas. Jitter estadístico, calibración absoluta óptica/material y reset eléctrico total son estudios posteriores, no requisitos pendientes del cierre aquí admitido.

## Reproducción y antecedentes

[Pipeline normal 05](../../../pipelines/05_effective_photon_run.py), [visual 05 de datos guardados](../../../plot_pipelines/05_plot_effective_photon_run.py) y [comandos presentes](../../GEMINGA_COMMANDS_PRESENT.md). Fuente única por ejecución, reanudación exacta, datos pesados/cachés en scratch y presupuesto global máximo del 90 %.

El [informe final](../../../output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf) y los tres documentos [Principal](../../../Matematica/Principal.pdf), [Auxiliar](../../../Matematica/Auxiliar.pdf), [Precomputo](../../../Matematica/Precomputo.pdf) separan resultados, ecuaciones principales 2D, auxiliares durante evolución y precomputados.

La [entrada anterior de etapa 5](../../releases/history/stage5_README_before_v1.1.0_20261012.md) conserva las decisiones de DC y del candidato Moyal rechazado. Su admisión DC no sustituye el nuevo DC efectivo; cada run05 prepara uno compatible con su energía. No se reabren campañas históricas por publicar este cierre. La etiqueta `v1.0.0` conserva la implementación de la memoria.
