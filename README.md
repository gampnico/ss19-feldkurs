# Processing Scintillometry Data in Complex Terrain

## 1. Scope of program

## 2. Feature documentation

### 2.1 Scintillometry Processing
### 2.2 Footprint Climatology 

## 3. Workflow

### 3.1 Scintillometry
### 3.2 Path Footprint Climatology

 - Create individual footprints for each point along the path length either
 by using the online 2D FFP tool, or the FFP_clim function provided by 
 Natascha Kljun.
 - Determine the resolution and cell size of each generated footprint via 
 the MATLAB engine for Python.
 - Create ASCII raster file (Python).
 - Generate TIFF files and apply weighting function to each TIFF.
 - Mosaic and average the TIFF files in R and generate the final contour plot.
 
 ### New procedure
 
 1. Generate footprints for entire path length using the FFP_clim function.
 2. Generate xllcenter, yllcenter coordinates for each footprint.
 2. Determine the resolution and cell size of each generated footprint.
 3. Calculate xllcorner,yllcorner coordinates:
 
    > `xllcorner = xllcenter - (nrow * cellsize)/2`
4. Generate ASCII raster files, inserting correct coordinates.
5. Generate TIFF files, apply weighting functions to each TIFF.
6. Mosaic and average TIFF files in R, generate final contour.
7. Layer contour plot over map (either QGIS or GIMP).