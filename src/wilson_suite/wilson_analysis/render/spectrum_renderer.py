"""
Use in WilsonSimulation:
self.rendering = renderer(self.spec, self.system, self.exp, self.diagn, self.name, self.spec_eval_setup)
self.rendering, self.diagn = renderer(self.spec, self.system, self.exp, self.diagn, self.name, self.spec_eval_setup)

Input:
self.spec, self.system, self.exp, self.diagn, self.name, self.spec_eval_setup; 

Steps: 
1. choice of projection - how to present spectrum data: 2d or 3d or 1d ... (figure dimension setup)
2. assignment of data to axes
3. dynamic range settings (for contour - colorbar, for 1d - threshold on y values?)
. style:
    a. axes titles
    b. axes ticks
    c. axes ticks labels
    d. axes labels
    e. title if any
. finalize figure styling: tight layout, etc...
. Save figure (where, filename)
"""
from abc import ABC, abstractmethod
from typing import Tuple, List, Any
import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ...wilson_main.spectrum_abstractions import EvaluationInfo, RenderingInfo


import logging
logger = logging.getLogger("wilson."+__name__)



class LevelCalculator:
    """
    Handles calculation of contour levels and normalization
    
    organizes the logic for computing levels and labels
    based on dynamic range, number of levels, and normalization type.
    """

    @staticmethod
    def compute_levels(intensities: float, dynamic_range: float,
                    nlevels: int, colormap_spacing: str = None,
                    reference_max: float = None) -> Tuple[np.ndarray, List[str]]:
        """Calculate levels for contours and colorbar ticks.
        
        If reference_max is provided, levels are computed relative to it
        instead of the data's own maximum.
        """
        if reference_max is not None:
            d_max = reference_max
            print(f'Scaling wrt reference_max: {reference_max:.2e}/{np.max(intensities):.2e}={reference_max/np.max(intensities):.2e}')
        else:
            d_max = np.max(intensities)
        
        if d_max <= 0:
            raise ValueError(
                "Logarithmic colormap requested, but data contains no positive values "
                f"(max={d_max})."
            )
        
        d_min = d_max / dynamic_range

        if colormap_spacing == "log":
            level_values = np.logspace(np.log10(d_min), np.log10(d_max), nlevels)
        
        elif colormap_spacing == "linear":
            level_values = np.linspace(d_min, d_max, nlevels)
        else:
            raise ValueError('Choose log or linear colormap_spacing')
        level_labels = [f"${val:.1e}$" for val in level_values]

        LevelCalculator._validate_levels(level_values)

        return level_values, level_labels

    @staticmethod
    def _validate_levels(levels):
        "copy from matplotlib"
        if len(levels) > 1 and np.min(np.diff(levels)) <= 0.0:
            raise ValueError("Contour levels must be increasing")
        
class SpectrumRenderer(ABC):
    """
    Abstract base class for spectrum rendering
    
    PlotConfig instance would normally be stored in the RenderingInfo instance
    but can be provided as config parameter

    spec_data - amplitudes data from evaluation procedure
    """
    
    def __init__(self, 
                #  result: EvaluatedResult, sealed: SealedSetup,
                 spec_data: np.ndarray | dict = None,
                 spec_grid: dict = None,
                 rnd_info: "RenderingInfo" = None, 
                 ev_info: "EvaluationInfo" = None,
                 do_diagn: bool = False):

        self.spec_data = spec_data
        self.rnd_info = rnd_info
        self.spec_grid = spec_grid

        # TODO not used currently
        self.do_diagn = do_diagn

        self.ev_info = ev_info

        # self.config = self.rnd_info.style_config
        self.level_calc = LevelCalculator()
        self.intensities = None

        self.Xdata = None
        self.Ydata = None
        self.Zdata = None


    @abstractmethod
    def initialize_plot(self) -> Any:
        """Initialize plotting surface"""
        pass

    @abstractmethod
    def create_contour(self, plot_obj: Any, levels: np.ndarray, data: np.ndarray) -> Any:
        """Create contour plot"""
        pass
    
    @abstractmethod
    def setup_axes(self, plot_obj: Any) -> Any:
        """Configure axes, ticks and labels"""
        pass

    @abstractmethod
    def add_colorbar(self, plot_obj: Any, levels: np.ndarray, labels: List[str]) -> None:
        """Add colorbar to plot"""
        pass
    
    @abstractmethod
    def finalize(self, plot_obj) -> None:
        """some finishing styling touches (positioning and resizing)"""
        pass

    @abstractmethod
    def save_plot(self, plot_obj: Any, filename: str) -> None:
        """Save plot to file"""
        pass
    
    def _create_data_masks(self, data: np.ndarray) -> Any:
        """
        """
        
        if self.rnd_info.apply_exp_magn_conditions_render:
            magn_conditions = self.rnd_info.exp_magn_conditions
        else:
            magn_conditions = None
        
        return compute_masks(data=data, 
                             dynamic_range=self.rnd_info.dynamic_range,
                             grid=self.spec_grid,
                             magn_conditions=magn_conditions, 
                             non_zero_margin=self.rnd_info.magn_conditions_margin)

    def prep_data(self, spec_data_operations: str) -> np.ndarray:
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

    def normalize_to_reference_max(self, reference_max: float=None):
        """
        not used right now
        """
        if reference_max is None:

            if self.rnd_info.reference_max is not None:
                reference_max = self.rnd_info.reference_max
            else:
                raise ValueError("Provide reference maximum value to normalize")

            self.intensities = self.intensities/reference_max

    def _validate_inputs(self):
        """
        data and settings should not contradict:
        - appropriate number of axes
        """
        if not isinstance(self.spec_data, np.ndarray):
            raise TypeError("spec_data should be a np.ndarray")
        if self.spec_data.size == 0:
            raise ValueError("spec_data array should not be empty")
        if not isinstance(self.spec_grid, dict):
            raise TypeError("spec_grid should be a dictionary with X,Y,(Z) data")

        for key, val in self.spec_grid.items():
            if not isinstance(val, np.ndarray):
                raise TypeError(f"spec_grid[{key!r}] is not a np.ndarray")
            if val.size == 0:
                raise ValueError(f"spec_grid[{key!r}] is an empty array")
        
        from ...wilson_main.spectrum_abstractions import EvaluationInfo, RenderingInfo
        if isinstance(self.rnd_info, RenderingInfo):
            self.config = self.rnd_info.style_config
        if not isinstance(self.ev_info, EvaluationInfo):
            raise TypeError("ev_info should be an instance of a class EvaluationInfo")
        if self.rnd_info.dynamic_range <= 0:
            raise ValueError("ev_info.dynamic_range must be positive")


    def _validate_data_2d(self):
        if self.Xdata is None:
            raise ValueError("Xdata was not set")
        if self.Ydata is None:
            raise ValueError("Ydata was not set")

        if self.intensities.ndim != 2:
            raise ValueError("intensities in renderer must be a 2D array")

        if self.Xdata.shape != self.intensities.shape or self.Ydata.shape != self.intensities.shape:
            raise ValueError("X,Y and intensities data do not match in shape:\n"
                             f"  x.shape = {self.Xdata.shape}\n"
                             f"  y.shape = {self.Ydata.shape}\n"
                             f"  z.shape = {self.intensities.shape}\n"
                             )

    def render(self, filename: str):
        """Main rendering pipeline"""
        self._validate_inputs()

        # prepare data for contour plotting with spec_data_operations and spec_grid.axes
        self.prep_data(spec_data_operations=self.rnd_info.spec_data_operations)
        self.normalize_to_reference_max(self.rnd_info.reference_max)
        self._validate_data_2d()

        # log10 = True if self.rnd_info.intensity_normalization_type is not None else False
        
        # Calculate levels with both original and normalized scales
        levels, labels = self.level_calc.compute_levels(
            intensities=self.intensities,
            dynamic_range=self.rnd_info.dynamic_range,
            nlevels=self.rnd_info.nlevels,
            colormap_spacing=self.config.colormap_spacing,
            reference_max=self.rnd_info.reference_max
        )
        
        self.levels = levels
        self.labels = labels

        fig, ax = self.initialize_plot()
        fig, ax, contour = self.create_contour(plot_obj=(fig, ax), levels=levels, data=self.intensities)
        fig, ax = self.setup_axes(plot_obj=(fig, ax))
        fig, ax, cbar = self.add_colorbar(plot_obj=(fig, ax, contour), levels=levels, labels=labels)
        
        self.finalize(plot_obj=(fig, ax, cbar))
        # self.save_plot(plot_obj=(fig, ax, cbar), filename=filename)

        return fig, ax, contour, cbar

def compute_masks(data: np.ndarray, 
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

