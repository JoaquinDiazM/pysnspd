# Referencia espacial térmica de un núcleo radial

Este ensayo prepara un contraste físico para D.4.4. El resultado buscado es la **derivada de energía libre respecto de la amplitud** sobre un mismo perfil prescrito. No busca una solución estacionaria del condensado ni una barrera de nucleación. No integra tiempo, no calibra KWT y no representa la geometría completa de la cinta de Korzh.

Se usa T=0,9 K, Tc=8,65 K y D=0,5 cm²/s. Las unidades son kBTc para amplitud y frecuencias de Matsubara, ℓ₀=√[ℏD/(2kBTc)]≈4,70 nm para distancias y N₀kBTc para la fuerza. La amplitud prescrita es d(r)=d_eq(T)tanh(r), con una vuelta de fase azimutal y potencial vector nulo. d_eq es la amplitud uniforme BCS a la temperatura impuesta.

## Problema espectral independiente

Para cada frecuencia εₙ=2π(T/Tc)(n+1/2), la simetría radial reduce Usadel a

```math
\theta_n''+\frac{\theta_n'}{r}
-\frac{\sin\theta_n\cos\theta_n}{r^2}
-\epsilon_n\sin\theta_n+d(r)\cos\theta_n=0.
```

La ecuación espectral de fase se satisface por simetría: la corriente azimutal no depende del ángulo. La rama buscada tiene 0≤θₙ<π/2 y es regular en el origen. Su expansión local es

```math
\theta_n(r)=c_n r+
\frac{\epsilon_n c_n-d'(0)-(2/3)c_n^3}{8}\,r^3+O(r^5).
```

En r_min=0,001 se emplea θ′ₙ=θₙ/r_min, el término principal de esa regularidad. El error de la condición es O(r_min²); no se declara una condición exacta en el origen ni se impone además θ′ₙ(0)=0. Esta última condición anularía una pendiente física que debe determinar el problema.

En r=R, θₙ se fija por el cierre uniforme **con circulación q=1/R**:

```math
\epsilon_n\sin\theta_n+R^{-2}\sin\theta_n\cos\theta_n
-d(R)\cos\theta_n=0.
```

Se resuelve su raíz física en tanθ: d=ε tanθ+Γ tanθ/√(1+tan²θ), una función estrictamente creciente. Este procedimiento conserva un intervalo de raíces y no recorta iteraciones de Newton para forzar una rama. La condición exterior sólo aproxima el espectro radial remoto; el contraste R=8 y R=12 mide su influencia en el mismo núcleo interior r≤4.

El radio R=12 corresponde a unos 56 nm, superior al semiancho de la cinta de 80 nm. Es admisible aquí porque se estudia un núcleo auxiliar circular de película extensa; no debe presentarse como dominio de simulación del dispositivo.

## Fuerza y cola espectral

La fuerza de referencia se calcula con

```math
X_U(r)=2d(r)\ln(T/T_c)
+4\pi(T/T_c)\sum_{n\ge0}\left[\frac{d(r)}{\epsilon_n}-\sin\theta_n(r)\right].
```

Su dependencia espacial ya está contenida en θₙ. **No se suma otro término −2K₀L₁d** a esta referencia: contaría de nuevo la respuesta de gradiente.

Las salidas principales conservan las sumas truncadas N=128 y N=256. También muestran, de manera separada, la estimación del primer término omitido. En el interior, definiendo L₁d=d″+d′/r−d/r²,

```math
\sin\theta_n=\frac d{\epsilon_n}+\frac{L_1d}{\epsilon_n^2}
+O(\epsilon_n^{-3}),\qquad
\Delta X_U^{\rm leading}=-4\pi(T/T_c)L_1d
\sum_{n\ge N}\epsilon_n^{-2}.
```

La suma omitida se evalúa con la función zeta de Hurwitz. Para tanh(r), L₁d es regular y comienza en −(8/3)d_eq r; el código usa una serie cerca de cero para evitar restar términos grandes. La estimación no se incorpora silenciosamente ni se trata como una cota rigurosa. Se comparan ambas versiones de las curvas, los dos cortes y los dos radios. La expansión interior no describe una capa artificial de borde; por eso los indicadores físicos se limitan a r≤4.

## Candidato comparado, en el mismo régimen estadístico

Se conserva la fuerza térmica del funcional regularizado C.9 calculada en C.7.2, para δ/Δ₀=0,05, 0,10 y 0,20. Su cierre uniforme incluye la resta analítica de la contribución lineal a gran Γ mediante digamma; el resto no lineal usa la cola de gran ε+Γ. Las raíces se calculan con bisección vectorial de una función monótona, sin el recorte de Newton del diagnóstico histórico.

Con K₀=π/4, ρ=d² y q_δ=ρ/[(ρ+δ²)r], la fuerza empleada es

```math
X_C=X_{\rm uni}-2K_0L_1d-2K_0dq_\delta^2
+(\Pi_{\rm uni}-2K_0\rho q_\delta)
\frac{2d\delta^2}{(\rho+\delta^2)^2r}.
```

δ se convierte a unidades kBTc mediante Δ₀/(kBTc)=πe^(−γ_E). El coeficiente de difusión distinto del diagnóstico histórico sólo cambia la conversión de r a nanómetros; la comparación adimensional conserva T/Tc y las mismas convenciones.

Ambos lados son derivados de **energía libre térmica a temperatura impuesta**. Los estados sintéticos a ocupación fija de los lotes 4A anteriores pertenecen a otro régimen estadístico y no se usan como referencia física de este ensayo. Se comparan perfiles, máximos y normas con peso 2πr. No se divide por fuerzas puntuales cercanas a cero. No se calcula energía total ni barrera: la cola logarítmica de energía del vórtice necesitaría un tratamiento común adicional.

## Ejecución y comprobaciones

El plan reutiliza los primeros 128 modos del caso R=8,N=256; hay 512 problemas espectrales distintos entre ambos radios y seis trabajos de cierre local. Todos comparten el mismo grupo de procesos, de modo que los casos liberan capacidad para los demás. El presupuesto de CPU y memoria se obtiene del mismo módulo que las campañas estáticas, con máximo del 90 %, afinidad registrada y una hebra BLAS por proceso. Hay barras, tiempo transcurrido y ETA. No hay reintentos, cambio silencioso de método ni sobrescritura de salidas.

Antes de la campaña se dispone de un piloto de sólo ocho modos. Las pruebas ligeras cubren una solución homogénea exacta sin circulación y una solución lineal independiente expresada mediante funciones de Bessel I₁/K₁, además de raíces con depareamiento fuerte y el límite regular del laplaciano complejo. El piloto no se interpreta como una fuerza física obtenida con ocho frecuencias.

Código: `sandbox/stage4_core/radial_usadel_reference.py`. Plan: `radial_reference_plan.json`. El piloto de ocho modos terminó en 1,19 s; con esa medida se admitió la campaña completa con un límite de 240 s. Terminó en 4,84 s, sin requerir entrega manual. El comando de reproducción se conserva en la libreta; la salida conserva identidad, campos espectrales, residuos del problema de borde, curvas y comparaciones.
