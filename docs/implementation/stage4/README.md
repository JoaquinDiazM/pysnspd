# Etapa 4 preparada: núcleo y disipación, primero sin fotón

El usuario autorizó cerrar la investigación 3.5 y preparar esta etapa con
estados controlados sin fotón. **La implementación de etapa 4 aún no comienza.**
El [contrato de entrada](entry_contract.json) y la
[revisión de capacidades reales](../stage3_5/assessment_r2_20260923/stage4_readiness.md)
definen el alcance. El [cierre de investigación 3.5](../stage3_5/assessment_r2_20260923/README.md)
conserva los datos, rangos condicionales y límites.

## Primera secuencia de trabajo

1. Exponer explícitamente los tiempos de KWT y el parámetro de núcleo en la
   ruta experimental, manteniendo sus valores heredados por defecto. Hoy los
   tiempos están escritos directamente en `KWTMobility`; editar un YAML no
   selecciona otro escenario. Guardar los valores realmente usados. Verificar
   las identidades de energía, fuerza, corriente y disipación que cambie esta ruta.
2. Preparar la referencia material coherente: Tc=8,65 K, Tb=0,9 K, D=0,5 cm²/s,
   R□=608 Ω y espesor7 nm; derivar N0 y conductividad conjuntamente. Conservar
   las tablas previas y distinguir sus funciones adimensionales de las escalas
   físicas del escenario nuevo. No volver a calcular tablas que puedan
   reutilizarse con una transformación de escalas justificada y registrada.
3. Registrar controles uniformes y perturbaciones suaves sin fotón, con
   poblaciones declaradas. Comparar energía, fuerza, corriente, signo de D.36
   y movilidad radial/tangencial. La comparación de δ/Δ0=0,05/0,10/0,20 prueba
   sensibilidad del cierre; no es una extrapolación automática a cero.
4. Cubrir cerca de cero la amplitud que necesite un núcleo. El catálogo positivo
   empieza en0,08Δ0 y su punto normal separado no permite interpolar ese hueco.
   Usar consulta causal directa o una extensión documentada, sin extrapolar
   la semilla rápida del estado casi uniforme.
5. Estimar el coste con un piloto acotado antes de preparar una campaña.
   Si el cálculo esperado supera cinco minutos, entregar el comando al usuario
   con progreso y ETA. No hay todavía una campaña larga lista para lanzar.

## Qué resultados se esperan

Un informe corto que distinga: identidades verificadas, sensibilidad a δ y
movilidad, frontera de estabilidad, soporte que falta y recomendación para
continuar o reformular. El error numérico debe permitir distinguir el cambio
físico estudiado; no se impone una tolerancia relativa universal sobre cantidades
casi nulas. Un símbolo principal negativo es una frontera del candidato, no
un error que pueda esconderse recortando rigidez o reduciendo el paso temporal.

## Después del diagnóstico inicial

Para dinámica espacial, completar el transiente débil y los balances de las
capacidades efectivamente usadas: poblaciones, potencial, circuito y bordes.
Una continuación1D requiere además su justificación transversal y el empalme
cinético; puede mantenerse el dominio enteramente2D mientras eso no exista.
El circuito conserva la topología y parámetros de la memoria, recalculando
la partición inductiva cuando cambien material o dominio.

El ancho fotónico, el reparto y reloj de transferencia, y las tasas absolutas
de NbN siguen sin admitirse para predecir hotbelt y latencia de Korzh. Esas
dependencias se resolverán antes de un ensayo que las use. No bloquean los
controles estáticos sin fotón. No se inicia jitter ni se ajusta el modelo para
forzar la diferencia experimental de4,2 ps.
