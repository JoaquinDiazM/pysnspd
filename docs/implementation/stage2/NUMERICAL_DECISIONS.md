# Registro de decisiones numéricas de la etapa 2

Los umbrales de `acceptance_criteria.json` se fijaron antes de los ensayos y se
mantienen. El anexo `review/complementary_grid_amendment.json` permite sustituir
la dimensión ocupacional de 180 estados por una malla complementaria, manteniendo
intactos el catálogo R2 y todos los umbrales. Los fallos no se reclasifican como
aceptados cuando se modifica la malla.

1. **Separar conservación de precisión.** La primera lista de pares conserva
   energía, pero una potencia de dispersión difiere un 0,1395 % y una respuesta
   de transporte un 0,2553 %, frente a 0,1 % permitido. Aumentar sólo el orden de
   integración no resuelve la representación de las poblaciones.
2. **Refinar las ocupaciones reales.** La malla complementaria utiliza el mismo
   potencial de vacío y la inversión causal directa del conteo. Tanto la fuerza
   como el trabajo proceden de esa energía. Cada malla inicializa de nuevo los
   perfiles analíticos declarados; no remapea silenciosamente un vector previo.
3. **Resolver las distribuciones además de sus potencias.** Una cuadratura en
   intervalos de conteo pasó las potencias pero falló la distribución fonónica.
   Alinear los intervalos con energías fonónicas corrigió ese fallo; todavía
   quedaba hasta 14,1 % en la respuesta electrónica. La partición conjunta debe
   resolver también los nodos electrónicos donde cambia el reparto del evento.
4. **Verificar la respuesta a la corriente.** Las energías y derivadas puntuales
   eran correctas, pero dos puntos de Gauss por intervalo cercano al borde
   dejaban un 1,56 % de error en un momento conjugado a Gamma. Refinar sólo la
   región de energías altas no lo corregía. La receta complementaria incorpora
   ocho puntos por intervalo cercano al borde y refina también esos intervalos.
5. **Usar la relación espectral exacta donde está disponible.** Para Gamma=0,
   los puntos virtuales de colisión pueden evaluarse con la inversión BCS causal
   cerrada. Se conserva el reparto positivo sobre los estados dinámicos y su
   energía; se evita aproximar innecesariamente la curvatura del espectro.

Los JSON de pilotos, sus configuraciones y los controles independientes se
conservan en esta carpeta. El dictamen final debe identificar el código, mallas y
trayectorias concretos que pasan, sin extender ese resultado a distribuciones
arbitrariamente estrechas, datos materiales sin calibrar o transientes espaciales.
