from typing import Tuple, List, Any
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from .spectrum_renderer import SpectrumRenderer 
from .render_utils import NormalizationType

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ...wilson_main.spectrum_abstractions import SpecEvalSetup
    from ...wilson_main.workflow_abstractions_updated import SealedSetup, EvaluatedResult

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pathlib import Path

import logging
logger = logging.getLogger("wilson."+__name__)

class MatplotlibRenderer(SpectrumRenderer):
    """Matplotlib implementation of spectrum renderer"""
    
    def initialize_plot(self) -> Tuple[plt.Figure, plt.Axes]:
        matplotlib.use('Agg')
        plt.rcParams['path.simplify'] = True
        plt.rcParams['agg.path.chunksize'] = 10000
        
        matplotlib.rc('font', **self.config.font_dict)
        
        fig = plt.figure(figsize=self.config.figsize)
        # Add axes with specific margins to ensure content fits
        ax = fig.add_axes([0.15, 0.15, 0.7, 0.75])  # [left, bottom, width, height]
        
        return fig, ax


    def create_contour(self, 
                       plot_obj: Tuple[plt.Figure, plt.Axes], 
                       levels: np.ndarray, 
                       data: np.ndarray) -> Tuple[plt.Figure, plt.Axes, Any]:
        """
        2D contour plot.

        data parameter here refers to the Z-axis data for contour plotting, the signal magnitude
        """

        fig, ax = plot_obj
        
        # Create masked arrays
        no_data_mask, below_range_mask = self._create_data_masks(data)

        if self.config.other_colors:
            # Fill below-range regions
            ax.contourf(self.Xdata, self.Ydata,
                    below_range_mask,
                    levels=[0, 0.5, 1],
                    colors=[self.config.below_range_color])

        # Setup base colormap
        cmap = plt.get_cmap(self.config.colormap).copy()
        cmap.set_over(self.config.saturation_color)
                
        if self.rnd_info.intensity_normalization_type is not None:
            # Create logarithmic normalization for color mapping
            norm = matplotlib.colors.LogNorm(vmin=levels[0], vmax=levels[-1])
        else:
            norm = None

        # Plot main data with normalized colors
        contour = ax.contourf(self.Xdata, self.Ydata, 
                           data,
                           levels=levels,
                           norm=norm,  # Add normalization
                           cmap=cmap,
                           extend='max')

        # Replace contourf with pcolormesh for the mask
        # Use a masked array so that 'False' values are completely transparent
        masked_no_data = np.ma.masked_where(~no_data_mask, no_data_mask)
        ax.contourf(self.Xdata, self.Ydata,
                masked_no_data,
                levels=[0, 0.5, 1],
                colors=[self.config.no_data_color])

        # Hide contour linestroke on pyplot.contourf to get only fills
        # https://stackoverflow.com/questions/8263769/hide-contour-linestroke-on-pyplot-contourf-to-get-only-fills
        contour.set_edgecolor("face")

        # if self.config.other_colors:   
        #     # Single clean edge line
        #     ax.contour(self.Xdata, self.Ydata,
        #             ~no_data_mask,
        #             levels=[0.5],
        #             colors=[self.config.data_edge_color],
        #             linewidths=[self.config.data_edge_width])

        return fig, ax, contour


    def setup_axes(self, plot_obj):
        fig, ax = plot_obj

        # Set up axes labels
        label_fontsize = self.config.label_fontsize if hasattr(self.config, 'label_fontsize') else 25

        xlabel_str = self.xyz_labels.get('x', r'$default x /2\pi c, \text{cm}^{-1}$')
        ylabel_str = self.xyz_labels.get('y', r'$default y /2\pi c, \text{cm}^{-1}$')

        # labelpad - distance from axis to label
        ax.set_xlabel(xlabel_str, fontsize=label_fontsize, labelpad=65.) 
        ax.set_ylabel(ylabel_str, fontsize=label_fontsize, labelpad=65.)


        # Simple aspect ratio setting
        if self.config.equal_aspect:
            ax.set_aspect('equal', adjustable='box')
        
        # Generate tick positions
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        
        # Set regular ticks
        ax.set_xticks(np.arange(
            np.ceil(x_min / self.config.tick_step) * self.config.tick_step,
            np.floor(x_max / self.config.tick_step) * self.config.tick_step + self.config.tick_step,
            self.config.tick_step
        ))
        
        ax.set_yticks(np.arange(
            np.ceil(y_min / self.config.tick_step) * self.config.tick_step,
            np.floor(y_max / self.config.tick_step) * self.config.tick_step + self.config.tick_step,
            self.config.tick_step
        ))
        
        # Set grid
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Set axis limits
        ax.set_xlim(
            self.config.x_min if self.config.x_min is not None else np.min(self.Xdata),
            self.config.x_max if self.config.x_max is not None else np.max(self.Xdata)
        )
        ax.set_ylim(
            self.config.y_min,
            self.config.y_max if self.config.y_max is not None else np.max(self.Ydata)
        )
        
        # After setting up the main axes ticks and labels
        if self.config.show_right_ticks:
            # Show ticks and labels on both sides of y-axis
            ax.yaxis.set_ticks_position('both')
            # Keep labels on both sides
            ax.tick_params(labelleft=True, labelright=True)
        
        if self.config.show_top_ticks:
            # Show ticks and labels on both sides of x-axis
            ax.xaxis.set_ticks_position('both')
            # Keep labels on both sides
            rotation=self.config.x_tick_rotation
            ax.tick_params(labelbottom=True, labeltop=True, axis='x', rotation=rotation)
        
        # Set axis limits
        ax.set_xlim(
            self.config.x_min if self.config.x_min is not None else np.min(self.Xdata),
            self.config.x_max if self.config.x_max is not None else np.max(self.Xdata)
        )
        ax.set_ylim(
            self.config.y_min,
            self.config.y_max if self.config.y_max is not None else np.max(self.Ydata)
        )
        
        # After setting ticks, rotate x-axis tick labels
        ax.tick_params(axis='x', rotation=self.config.x_tick_rotation)
        # https://stackoverflow.com/questions/2969867/how-do-i-add-space-between-the-ticklabels-and-the-axes
        
        return fig, ax
    
    def add_colorbar(self, plot_obj: Tuple[plt.Figure, plt.Axes, Any], 
                    levels: np.ndarray, labels: List[str]) -> Tuple[plt.Figure, plt.Axes, Any]:
        """
        https://pythonmatplotlibtips.blogspot.com/2019/07/draw-two-axis-to-one-colorbar.html
        """
        fig, ax, contour = plot_obj

        # Create colorbar and manually align it to the plot's height
        cbar = fig.colorbar(contour, ax=ax)

        ax1 = cbar.ax
        ax1.set_aspect('auto')

        fig.canvas.draw()  # Ensure layout is updated
        pos = cbar.ax.get_position()

        # Create and set up normalized (left) axis
        ax2 = ax1.twinx()
        ax2.set_position(pos)
        
        # Calculate normalized positions based on selected normalization type
        if self.rnd_info.intensity_normalization_type == NormalizationType.LOG_RATIO:
            norm_positions = np.log10(levels)/np.log10(levels[-1])
            norm_format = "{x:.3f}"
            norm_label = "Log Ratio"
        elif self.rnd_info.intensity_normalization_type == NormalizationType.DECIBEL:
            norm_positions = 10 * np.log10(levels/levels[-1])
            norm_format = "{x:.1f} dB"
            norm_label = "Intensity (dB)"
        elif self.rnd_info.intensity_normalization_type == NormalizationType.PERCENTAGE:
            norm_positions = (levels/levels[-1]) * 100
            norm_format = "{x:.1f}%"
            norm_label = "Relative Intensity (%)"
        elif self.rnd_info.intensity_normalization_type is None:
            norm_positions = levels
            norm_format = "{x:.2f}"
            norm_label = "Original"

        else:  # LOG_SCALE
            norm_positions = (np.log10(levels) - np.log10(levels[0]))/(np.log10(levels[-1]) - np.log10(levels[0]))
            norm_format = "{x:.2f}"
            norm_label = "Log-scale Normalized"
        
        if self.rnd_info.intensity_normalization_type is not None:
            logger.debug(f"Normalized positions ({self.rnd_info.intensity_normalization_type.value}): {norm_positions}") #z
        
        # Set up normalized axis limits and ticks
        ax2.set_ylim(min(norm_positions), max(norm_positions))
        ax2.set_yticks(norm_positions)
        ax2.set_yticklabels([norm_format.format(x=x) for x in norm_positions])
        ax2.yaxis.set_ticks_position('left')
        ax2.yaxis.set_label_position('left')
        
        # Move colorbar position slightly to the right
        pos.x0 += 0.06
        pos.x1 += 0.06
        
        # Set up main (right) axis
        ax1.set_position(pos)
        ax1.yaxis.set_ticks_position('right')
        ax1.yaxis.set_label_position('right')
        ax1.set_yticks(levels)
        ax1.set_yticklabels(labels)
        ax1.set_ylabel(self.config.colorbar_main_label,
                       rotation=90,
                       labelpad=48, # distance from axis to label
                       fontsize=self.config.label_fontsize if hasattr(self.config, 'label_fontsize') else 25)
        
        # Set normalized axis label
        ax2.set_ylabel(norm_label,
                       rotation=90,
                       labelpad=48, # distance from axis to label
                       fontsize=self.config.label_fontsize if hasattr(self.config, 'label_fontsize') else 25)
        
        # Adjust spacing between axes
        cbar.ax.spines['right'].set_position(('outward', 0))
        ax2.spines['left'].set_position(('outward', 0))
        return fig, ax, cbar

    def finalize(self, plot_obj: Tuple[plt.Figure, plt.Axes, Any]) -> None:
        """
        Finalize plot styling.
        This method can be overridden in subclasses for additional styling.
        """
        
        fig, ax, cbar = plot_obj
        fig.canvas.draw()  # Make sure the figure layout is updated

        # Get the main axis position (where the actual plot is)
        ax_pos = ax.get_position()

        if self.rnd_info.style_config.axes_limits is not None:
            ax.set_xlim(*self.rnd_info.style_config.axes_limits['x'])
            ax.set_ylim(*self.rnd_info.style_config.axes_limits['y'])

        # Example: shrink or stretch colorbar to match ax height
        cbar_pos = cbar.ax.get_position()
        cbar.ax.set_position([
            cbar_pos.x0+self.config.colorbar_padding,       # x-position (keep same or adjust)
            ax_pos.y0,         # align bottom of colorbar to ax
            cbar_pos.width,    # keep same width
            ax_pos.height      # match ax height
        ])

    def save_plot(self, plot_obj: Tuple[plt.Figure, plt.Axes, Any], filename: str) -> None:
        fig, ax, _ = plot_obj
        
        # No need to set aspect here anymore as it's handled in create_contour
        ax.grid(True, linestyle='--', alpha=0.7)
        fig.savefig(filename, bbox_inches='tight',
                    dpi=self.config.dpi, format=filename.split('.')[-1])
        plt.close(fig)

def spectral_axis_to_label(axis_dict: dict, divide_by_2pic: bool = True) -> str:
    """
    Utility function.
    Making labels for axes using SpectralAxis.freq_vars dict
    
    axis_dict = SpectralAxis.freq_vars
    """
    terms = []
    for key, coeff in sorted(axis_dict.items(), key=lambda x: x[0]):
        if coeff == 0:
            continue
        sign = '+' if coeff > 0 else '-'
        abs_coeff = abs(coeff)

        if abs_coeff == 1:
            term = f"{sign} \\omega_{{{key[1:]}}}"  # remove 'w' from 'w1'
        else:
            term = f"{sign} {abs_coeff}\\omega_{{{key[1:]}}}"
        terms.append(term)

    if not terms:
        expr = "0"
    else:
        expr = " ".join(terms)
        if expr.startswith('+ '):
            expr = expr[2:]  # remove leading '+ ' for aesthetics

    # Wrap in parentheses if there are multiple terms
    if len(terms) > 1:
        expr = f"({expr})"

    if divide_by_2pic:
        expr = rf"${expr}/2\pi c, \text{{cm}}^{{-1}}$"
    else:
        expr = rf"${expr}, \text{{cm}}^{{-1}}$"

    return expr



class MatplotlibRendererGV:
    """
    
    """
    def __init__(self, 
                 eval_result: 'EvaluatedResult', 
                 setup_inputs: 'SealedSetup'):
        
        self.eval_result = eval_result
        self.setup_inputs = setup_inputs

        self.config = self.setup_inputs.spec_eval.rnd_info.style_config

    # -- conventional view methods --

    def _init_plot(self) -> Tuple[plt.Figure, plt.Axes]:
        matplotlib.use('Agg')
        plt.rcParams['path.simplify'] = True
        plt.rcParams['agg.path.chunksize'] = 10000
        
        matplotlib.rc('font', **self.config.font_dict)
        
        fig = plt.figure(figsize=self.config.figsize)
        # Add axes with specific margins to ensure content fits
        ax = fig.add_axes([0.15, 0.15, 0.7, 0.75])  # [left, bottom, width, height]
        
        return fig, ax

    def _prep_data(self, spec_data_operations: str) -> np.ndarray:
        """
        Prepare data for contour plotting.
        """
        if self.intensities is None:
            if spec_data_operations == 'abs()**2':
                self.intensities = np.abs(self.spec_data) ** 2
            elif spec_data_operations == 'abs':
                self.intensities = np.abs(self.spec_data)
            elif spec_data_operations == 'real':
                self.intensities = np.real(self.spec_data)
            elif spec_data_operations == 'imag':
                self.intensities = np.imag(self.spec_data)
            elif spec_data_operations == 'none':
                self.intensities = self.spec_data
            else:
                raise ValueError(f"Unsupported spec_data_operations: {spec_data_operations}")
        
        if self.spec_grid is None:
            raise ValueError('This SpectrumRenderer.spec_grid is None')        

        self.xyz_labels = {'x': None, 'y': None, 'z': None}
        if self.rnd_info.axes_labels is None:
            for i, o_k in enumerate(list(self.spec_grid.keys())):
                self.xyz_labels[list(self.xyz_labels.keys())[i]] = o_k
        else:
            self.xyz_labels = self.rnd_info.axes_labels
        
        if len(self.spec_grid)==3:
            self.Xdata, self.Ydata, self.Zdata = list(self.spec_grid.values())
        elif len(self.spec_grid)==2:
            self.Xdata, self.Ydata = list(self.spec_grid.values())

    def _normalize_to_reference_max(self, reference_max: float=None):
        """
        not used right now
        """
        if reference_max is None:

            if self.rnd_info.reference_max is not None:
                reference_max = self.rnd_info.reference_max
            else:
                raise ValueError("Provide reference maximum value to normalize")

            self.intensities = self.intensities/reference_max

    def _compute_masks(self, data: np.ndarray, 
                    dynamic_range: float, 
                    grid: dict=None,
                    magn_conditions: tuple[tuple]=None,
                    non_zero_margin: float=80.):
        """
        Intensities should be > 0, not negative

        data: intensities array
        non_zero_margin: added to 0, which is the boundary in magn_conditions
        """
        # ONLY EVV w2>w1 for paper 1 now
        if magn_conditions==(('B',),):
            if grid is None:
                raise ValueError("in compute_masks() grid is None")
            no_data = grid['B'] < (0 + non_zero_margin)
        elif magn_conditions==(('-A', 'B',),):
            if grid is None:
                raise ValueError("in compute_masks() grid is None")
            no_data = grid['B'] - grid['A'] < (0 + non_zero_margin)
        else:
            no_data = np.isnan(data)

        d_max = np.nanmax(data)

        if not np.isfinite(d_max) or d_max <= 0:
            below = np.zeros_like(data, dtype=bool)
        else:
            d_min = d_max / dynamic_range
            below = (~no_data) & (data < d_min)

        return no_data, below

    def _create_data_masks(self, data: np.ndarray) -> Any:
        """
        """
        
        if self.rnd_info.apply_exp_magn_conditions_render:
            magn_conditions = self.rnd_info.exp_magn_conditions
        else:
            magn_conditions = None
        
        return self._compute_masks(data=data, 
                                   dynamic_range=self.rnd_info.dynamic_range,
                                   grid=self.spec_grid,
                                   magn_conditions=magn_conditions, 
                                   non_zero_margin=self.rnd_info.magn_conditions_margin)


    # def contour(self, *, slice_spec=None, canvas=None, ) -> 'Axes':
    def contour(self, *, 
                slice_spec: dict | None = None, 
                ax: Axes | None = None, **opts) -> 'Axes':

        """
        2D contour plot. Overlays on canvas if given.
        If spec_data is >2D, pass slice_spec (e.g., {'dim2': 0.5}) to reduce.
        
        opts:
            cmap: str = 'viridis', levels: int = 20,

        1. prep Z data
        2. set up fig and ax?
        3. make levels
        4. cmap in contourf
        5. norm in contourf
        6. make colorbar
        7. ticks and axes customization
        """
        # 1. prep Z data
        data_2d = self._reduce_to_2d(slice_spec)
        data_2d = np.abs(data_2d)**2

        # 2. set up fig and ax?
        if ax is None:
            fig, ax = plt.subplots(figsize=self.config.figsize)
        else:
            fig = ax.figure

        # 3. make levels
        from .spectrum_renderer import LevelCalculator
        levels, labels = LevelCalculator().compute_levels(
            intensities=data_2d,
            dynamic_range=self.setup_inputs.spec_eval.rnd_info.dynamic_range,
            nlevels=self.setup_inputs.spec_eval.rnd_info.nlevels,
            colormap_spacing=self.config.colormap_spacing,
            reference_max=self.setup_inputs.spec_eval.rnd_info.reference_max
        )
        
        # 4. cmap in contourf
        cmap = plt.get_cmap(self.config.colormap).copy()
        cmap.set_over(self.config.saturation_color)

        # 5. norm in contourf
        if self.setup_inputs.spec_eval.rnd_info.intensity_normalization_type is not None:
            # Create logarithmic normalization for color mapping
            norm = matplotlib.colors.LogNorm(vmin=levels[0], vmax=levels[-1])
        else:
            norm = None

        # Plot main data with normalized colors
        contour = ax.contourf(self.eval_result.grid['A'], 
                              self.eval_result.grid['B'], 
                              data_2d,
                              levels=levels,
                              norm=norm,  # Add normalization
                              cmap=cmap,
                              extend='max')

        # Hide contour linestroke on pyplot.contourf to get only fills
        # https://stackoverflow.com/questions/8263769/hide-contour-linestroke-on-pyplot-contourf-to-get-only-fills
        contour.set_edgecolor("face")

        # 6. make colorbar
        fig, ax, cbar = add_colorbar((fig, ax, contour), 
                                     self.setup_inputs.spec_eval.rnd_info,
                                     self.config, levels, labels)
        ax_tick_fs = getattr(self.config, 'tick_label_fontsize', 25)
        ax.tick_params(axis='x', labelsize=ax_tick_fs)
        ax.tick_params(axis='y', labelsize=ax_tick_fs)

        return ax


    def line(self, *, dim, at, canvas=None, **opts) -> 'Axes':
        """1D line (slice of spec_data along `dim` at `at`)."""
        raise NotImplementedError

    def scatter_features(self, *, features=None, label=False, canvas=None) -> 'Axes':
        """Feature markers. Defaults to self._result.features."""
        raise NotImplementedError

    def scatter_plt(self, *, points=None, label=False, canvas=None) -> 'Axes':
        """Feature markers. Defaults to self._result.features."""
        raise NotImplementedError

    # -- conventional composite --

    def contour_with_features(self, **opts) -> 'Figure':
        """Contour + feature overlay on one figure (common case)."""
        raise NotImplementedError


    def save(self, obj, path: Path | str, **kwargs) -> Path:
        """Save a returned object to disk. Accepts Figure, Axes, DataFrame."""
        if hasattr(obj, 'savefig'):              # matplotlib Figure
            obj.savefig(path, **kwargs)
        elif hasattr(obj, 'get_figure'):          # matplotlib Axes
            obj.get_figure().savefig(path, **kwargs)
        elif hasattr(obj, 'to_csv'):              # DataFrame
            obj.to_csv(path, **kwargs)
        else:
            raise TypeError(f"Don't know how to save {type(obj).__name__}")
        return Path(path)

    def _reduce_to_2d(self, slice_spec):
        """Handle N-D → 2D reduction. slice_spec maps dim name to value."""
        return self.eval_result.spec


# def add_colorbar(plot_obj: Tuple[plt.Figure, plt.Axes, Any],
#                  rnd_info, config,
#                  levels: np.ndarray, labels: List[str]) -> Tuple[plt.Figure, plt.Axes, Any]:
#     """
#     https://pythonmatplotlibtips.blogspot.com/2019/07/draw-two-axis-to-one-colorbar.html
#     """
#     fig, ax, contour = plot_obj

#     # Create colorbar and manually align it to the plot's height
#     cbar = fig.colorbar(contour, ax=ax)

#     ax1 = cbar.ax
#     ax1.set_aspect('auto')

#     fig.canvas.draw()  # Ensure layout is updated
#     pos = cbar.ax.get_position()

#     # Create and set up normalized (left) axis
#     ax2 = ax1.twinx()
    
#     # Calculate normalized positions based on selected normalization type
#     if rnd_info.intensity_normalization_type == NormalizationType.LOG_RATIO:
#         norm_positions = np.log10(levels)/np.log10(levels[-1])
#         norm_format = "{x:.3f}"
#         norm_label = "Log Ratio"
#     elif rnd_info.intensity_normalization_type == NormalizationType.DECIBEL:
#         norm_positions = 10 * np.log10(levels/levels[-1])
#         norm_format = "{x:.1f} dB"
#         norm_label = "Intensity (dB)"
#     elif rnd_info.intensity_normalization_type == NormalizationType.PERCENTAGE:
#         norm_positions = (levels/levels[-1]) * 100
#         norm_format = "{x:.1f}%"
#         norm_label = "Relative Intensity (%)"
#     elif rnd_info.intensity_normalization_type is None:
#         norm_positions = levels
#         norm_format = "{x:.2f}"
#         norm_label = "Original"

#     else:  # LOG_SCALE
#         norm_positions = (np.log10(levels) - np.log10(levels[0]))/(np.log10(levels[-1]) - np.log10(levels[0]))
#         norm_format = "{x:.2f}"
#         norm_label = "Log-scale Normalized"
    
#     if rnd_info.intensity_normalization_type is not None:
#         logger.debug(f"Normalized positions ({rnd_info.intensity_normalization_type.value}): {norm_positions}") #z
    
#     # Set up normalized axis limits and ticks
#     ax2.set_ylim(min(norm_positions), max(norm_positions))
#     ax2.set_yticks(norm_positions)
#     ax2.set_yticklabels([norm_format.format(x=x) for x in norm_positions])
#     ax2.yaxis.set_ticks_position('left')
#     ax2.yaxis.set_label_position('left')
    
#     # Move colorbar position slightly to the right
#     pos.x0 += 0.06
#     pos.x1 += 0.06
    
#     # Set up main (right) axis
#     ax2.set_position(pos)
#     ax1.set_position(pos)
#     ax1.yaxis.set_ticks_position('right')
#     ax1.yaxis.set_label_position('right')
#     ax1.set_yticks(levels)
#     ax1.set_yticklabels(labels)
#     ax1.set_ylabel(config.colorbar_main_label,
#                     rotation=90,
#                     labelpad=48, # distance from axis to label
#                     fontsize=config.label_fontsize if hasattr(config, 'label_fontsize') else 25)
    
#     # Set normalized axis label
#     ax2.set_ylabel(norm_label,
#                     rotation=90,
#                     labelpad=48, # distance from axis to label
#                     fontsize=config.label_fontsize if hasattr(config, 'label_fontsize') else 25)
    
#     # Adjust spacing between axes
#     cbar.ax.spines['right'].set_position(('outward', 0))
#     ax2.spines['left'].set_position(('outward', 0))


#     ax_pos = ax.get_position()
#     if rnd_info.style_config.axes_limits is not None:
#         ax.set_xlim(*rnd_info.style_config.axes_limits['x'])
#         ax.set_ylim(*rnd_info.style_config.axes_limits['y'])

#     cbar_pos = cbar.ax.get_position()

#     cbar.ax.set_position([
#         cbar_pos.x0+config.colorbar_padding,       # x-position (keep same or adjust)
#         ax_pos.y0,         # align bottom of colorbar to ax
#         cbar_pos.width,    # keep same width
#         ax_pos.height      # match ax height
#     ])


#     final_rect = [
#         cbar_pos.x0 + config.colorbar_padding,
#         ax_pos.y0,
#         cbar_pos.width,
#         ax_pos.height,
#     ]
#     cbar.ax.set_position(final_rect)
#     ax2.set_position(final_rect)
    
#     return fig, ax, cbar



from matplotlib.ticker import FixedLocator, FixedFormatter

def add_colorbar(plot_obj: Tuple[plt.Figure, plt.Axes, Any],
                 rnd_info, config,
                 levels: np.ndarray, labels: List[str]) -> Tuple[plt.Figure, plt.Axes, Any]:
    """
    Colorbar with real level values on the right and normalized values on the
    left. Left axis is locked to the right via secondary_yaxis, so alignment
    is guaranteed regardless of the cbar's internal norm.
    """
    fig, ax, contour = plot_obj
    assert np.all(np.diff(levels) > 0), "levels must be strictly increasing"

    ref = levels[-1]
    log_lo, log_hi = np.log10(levels[0]), np.log10(ref)
    EPS = 1e-300

    # (forward, inverse, format, label) per normalization type
    NT = NormalizationType
    norm_specs = {
        NT.LOG_RATIO:  (lambda y: np.log10(np.maximum(y, EPS)) / log_hi,
                        lambda y: 10.0 ** (y * log_hi),
                        "{x:.3f}",   "Log Ratio"),
        NT.DECIBEL:    (lambda y: 10 * np.log10(np.maximum(y, EPS) / ref),
                        lambda y: ref * 10.0 ** (y / 10),
                        "{x:.1f} dB", "Intensity (dB)"),
        NT.PERCENTAGE: (lambda y: 100 * y / ref,
                        lambda y: y * ref / 100,
                        "{x:.1f}%",  "Relative Intensity (%)"),
        NT.LOG_SCALE:  (lambda y: (np.log10(np.maximum(y, EPS)) - log_lo) / (log_hi - log_lo),
                        lambda y: 10.0 ** (log_lo + y * (log_hi - log_lo)),
                        "{x:.2f}",   "Log-scale Normalized"),
        None:          (lambda y: y, lambda y: y, "{x:.2f}", "Original"),
    }
    forward, inverse, norm_fmt, norm_label = norm_specs[rnd_info.intensity_normalization_type]
    norm_positions = forward(levels)

    # --- Colorbar ---
    cbar = fig.colorbar(contour, ax=ax)
    ax1 = cbar.ax
    ax1.set_aspect('auto')

    if rnd_info.style_config.axes_limits is not None:
        ax.set_xlim(*rnd_info.style_config.axes_limits['x'])
        ax.set_ylim(*rnd_info.style_config.axes_limits['y'])
    fig.canvas.draw()

    # Pin cbar to main-axes height
    ax_pos, cbar_pos = ax.get_position(), ax1.get_position()
    ax1.set_position([cbar_pos.x0 + 0.06 + config.colorbar_padding,
                      ax_pos.y0, cbar_pos.width, ax_pos.height])

    label_fs = getattr(config, 'label_fontsize', 25)
    tick_fs  = getattr(config, 'cb_tick_label_fontsize', 14)

    # Right axis: real level values
    cbar.set_ticks(levels)
    cbar.set_ticklabels(labels)
    ax1.tick_params(axis='y', labelsize=tick_fs)
    ax1.set_ylabel(config.colorbar_main_label, rotation=90, labelpad=48, fontsize=label_fs)

    # Left axis: normalized values, auto-aligned
    ax2 = ax1.secondary_yaxis('left', functions=(forward, inverse))
    ax2.yaxis.set_major_locator(FixedLocator(norm_positions))
    ax2.yaxis.set_major_formatter(
        FixedFormatter([norm_fmt.format(x=x) for x in norm_positions])
    )
    ax2.tick_params(axis='y', labelsize=tick_fs)
    ax2.set_ylabel(norm_label, rotation=90, labelpad=48, fontsize=label_fs)

    return fig, ax, cbar

