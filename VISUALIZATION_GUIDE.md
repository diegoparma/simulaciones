# 3D Spin Visualization Tools - Complete Guide

## Overview

Te proporcioné tres herramientas para visualizar la dinámica de spines en 3D:

### 1. **Plotly/HTML (Recomendado - Más rápido)**
- **Archivo**: `visualize_spins_plotly.py`
- **Pros**: 
  - Interactivo, rápido, en el navegador
  - No requiere instalación de OVITO
  - Genera archivos HTML portátiles
  - Zoom, rotación, hover tooltip con valores de spin
- **Cons**: Menos features avanzadas que OVITO

### 2. **OVITO GUI (Profesional)**
- **Archivo**: `open_in_ovito.py`
- **Pros**: 
  - Herramienta de visualización científica profesional
  - Vector display con propiedades derivadas
  - Análisis avanzado (estructura, correlaciones)
  - Publicable
- **Cons**: Interfaz más compleja, requiere instalación

### 3. **OVITO Script (Programático)**
- **Archivo**: `ovito_visualize_spins.py`
- **Pros**: 
  - Automatización completa
  - Batch processing
  - Exportación a múltiples formatos
- **Cons**: Requiere OVITO Pro para rendering programático

---

## Quick Start

### Opción A: Visualización rápida en HTML (RECOMENDADO)

```bash
cd ~/Documents/simulaciones
source venv_viz/bin/activate
python3 visualize_spins_plotly.py
```

Se crearán archivos (abre cualquiera en tu navegador):
- `spin_frame_000.html` - Primer frame (equilibrio)
- `spin_frame_025.html` - Frame intermedio (SOT activo)
- `spin_frame_050.html` - Último frame (máximo efecto)
- `spin_animation_latest.html` - Vista final en detalle

En el navegador:
- **Click + Drag** = Rotar vista
- **Scroll** = Zoom in/out
- **Hover** = Ver spin vector de cada átomo

---

### Opción B: Abrir en OVITO GUI

```bash
cd ~/Documents/simulaciones
python3 open_in_ovito.py
```

O directamente:
```bash
/Applications/Ovito.app/Contents/MacOS/ovito dump.spins.lammpstrj
```

En la ventana de OVITO:
1. Verás la estructura atómica en 3D
2. Click en diferentes frames en la timeline (abajo)
3. Para ver vectores spin:
   - Add Modifier → Vector (busca)
   - Configure para mostrar c_outsp[1,2,3]
   - Aumenta tamaño/color para visibilidad

---

### Opción C: Script OVITO programático (Avanzado)

```bash
cd ~/Documents/simulaciones
ovitos ovito_visualize_spins.py
```

Exporta varias visualizaciones y formatos automáticamente.

---

## Interpretación de las visualizaciones

### Colores en los átomos
- **Brillo/Color**: Magnitud del spin
- **Rojo**: Spins débiles
- **Azul/Verde**: Spins fuertes
- Hover muestra: `Atom X: S=(sx, sy, sz)`

### Vectores de spin (flechas)
- **Largo**: Magnitud del spin
- **Orientación**: Dirección de magnetización
- **Comparar frames**:
  - Frame 0 (equilibrio): Spins aleatorios/débiles
  - Frame 25-50 (con SOT): Spins orientados, cambian con tiempo

---

## Ejemplos de análisis

### Ver efecto SOT en tiempo real

```bash
# Frame inicial (antes de SOT)
open spin_frame_000.html

# Frame SOT activado
open spin_frame_025.html

# Cambios evidentes:
# - Sy aumentó (vectores giran en plano xy)
# - |S| creció (spines se alinearon más)
```

### Comparar configuraciones

```bash
# Lado a lado en dos ventanas del navegador:
open spin_frame_000.html &
open spin_frame_050.html &
```

Compara visual:
- Distribución de vectores
- Orientación preferida
- Correlaciones espaciales

---

## Para documentar/presentar resultados

### Opción 1: Screenshots desde HTML
1. Abre en Chrome/Safari
2. Click derecho → "Take Screenshot"
3. Captura para presentaciones

### Opción 2: Exportar desde OVITO
1. En OVITO: File → Export Video
2. Selecciona rango de frames
3. Genera MP4/MOV de la simulación

### Opción 3: Datos numéricos
Ya tenemos: `analyze_spins_simple.py`
```bash
source venv_viz/bin/activate
python3 analyze_spins_simple.py dump.spins.lammpstrj
```

Produce tabla con cambio de magnetización (ya ejecutado, mostró +0.55% en Sy).

---

## Troubleshooting

### "ModuleNotFoundError: plotly"
```bash
cd ~/Documents/simulaciones
source venv_viz/bin/activate  # Activa el virtual env
pip install plotly
```

### "OVITO not found"
```bash
brew install ovito
```

### Archivos .html vacíos o lentos
- Reduce el número de frames (modifica `visualize_spins_plotly.py`)
- O usa versión simplificada con fewer particles

### No veo los vectores de spin en OVITO

Asegúrate:
1. Pipeline tiene DATA SOURCE con LAMMPS dump
2. Agregaste Vector Display modifier
3. En Vector Display settings, mapeaste las propiedades:
   - Vector X: c_outsp[1]
   - Vector Y: c_outsp[2]
   - Vector Z: c_outsp[3]

---

## Flujo recomendado

```
1. visualize_spins_plotly.py 
   ↓
2. Revisar spin_frame_*.html en navegador
   ↓
3. Si todo OK, usar open_in_ovito.py para análisis detallado
   ↓
4. Exportar figures/videos desde OVITO para paper/presentación
```

---

## Archivos generados

```
After running visualize_spins_plotly.py:

spin_frame_000.html          (51 KB) - Inicio
spin_frame_025.html          (51 KB) - Medio  
spin_frame_050.html          (51 KB) - Final
spin_animation_latest.html   (51 KB) - Detalle último frame
venv_viz/                           - Virtual environment (3.2 GB)
```

Archivos seguros para compartir: todos los `.html` 
- Son independientes (contienen datos + javascript)
- Funcionan offline
- Abren en cualquier navegador

---

## Próximos pasos

Con esta visualización podés:

✅ **Verificar SOT funciona**: Ver cambios en spin orientación
✅ **Publicar resultados**: Figuras HTMLs o videos OVITO
✅ **Debugging**: Verificar si hay artefactos numéricos
✅ **Presentación**: Screenshots de física de spines

Luego podrías:
- 🎯 Hacer prueba de **skyrmión real** (DMI + SOT)
- 📊 Crear diagrama **velocidad vs corriente**
- 💻 Optimizar parámetros para **máximo drift**
