# wilson-analysis


## Rendering

https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html


### Assumptions about the spectrum intensity values:
- order of magnitude 10**5..8 ---> log spaced contour levels for contour 2D plots


### Rendering workflow

1. Prepare data
2. Plot
    - Contour plot - a 2D plot with Z (intensity) data as a 2D array.
  
        Extra prep:
        - levels
        - normalization
        - cmap choice/customization
        - colorbar setup
    - Line plot with Y (intensity)
    - Scatter plot with X and Y as position (X,Y), possibly with color or marker size as other (intensity) dimensions
    - (Analyses) Statistical plots: bar charts, histograms
3. Customizations:
    - axes
    - colors
    - markers
    - ticks
    - labels
    - title
    - text labels/annotations

