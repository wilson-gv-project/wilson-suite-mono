from enum import Enum
from typing import Optional, Tuple, Dict, Any
from dataclasses import dataclass, field


class NormalizationType(Enum):
    """
    # LOG_RATIO: log10(x)/log10(max)
    # Shows relative order of magnitude
    # Example: 0.5 means halfway between min and max in log scale

    # DECIBEL: 10 * log10(x/max)
    # Standard intensity scale in spectroscopy
    # 0 dB = max, -20 dB = 1/100 of max
    https://www.montana.edu/rmaher/eele417_fl14/decibel_scale_eele417.pdf
    https://en.wikipedia.org/wiki/Decibel
    https://www.animations.physics.unsw.edu.au/jw/dB.htm
    
    # PERCENTAGE: (x/max) * 100
    # Linear scale percentage
    # Direct proportion to maximum

    # LOG_SCALE: (log10(x) - log10(min))/(log10(max) - log10(min))
    # Normalized position in log space
    # 0 = minimum, 1 = maximum
    """
    LOG_RATIO = "log_ratio"
    DECIBEL = "db"
    PERCENTAGE = "percent"
    LOG_SCALE = "log_scale"

@dataclass
class PlotConfig:
    """Configuration for plot styling"""
    figsize: Tuple[int, int] = (35, 45)
    label_fontsize: int = 25
    cb_tick_label_fontsize: int = 15
    tick_label_fontsize: int = 15
    font_dict: Dict[str, Any] = field(default_factory=lambda: {'size': 20})
    colormap: str = 'magma'  # Better contrast colormap
    saturation_color: str = '#FF00FF'
    dpi: int = 250
    tick_step: float = 200.0  # Step size for both axes ticks
    equal_aspect: bool = True  # Force equal aspect ratio for axes
    other_colors: bool = True
    no_data_color: str = '#E0E0E0'  # Light gray
    below_range_color: str = '#F8F8F8'  # Very light gray
    data_edge_color: str = 'black'
    data_edge_width: float = 0.75
    x_min: Optional[float] = None
    x_max: Optional[float] = None
    y_min: Optional[float] = None
    y_max: Optional[float] = None
    colorbar_main_label: str = "Intensity"
    colorbar_padding: float = 0.01  # Padding between colorbar and plot
    show_top_ticks: bool = False
    show_right_ticks: bool = False
    x_tick_rotation: float = 45  # Add this line for configurable rotation
    colormap_spacing: str = "log"  # Options: "log", "linear"
    axes_limits: dict = None

    def __post_init__(self):
        if not isinstance(self.figsize, tuple):
            raise TypeError("figsize needs to be given as a tuple")
        if any(i < 0 for i in self.figsize):
            raise ValueError("Negative figsize was provided")
        if self.tick_step < 0:
            raise ValueError("Negative tick_step was provided")
        
        if self.x_min is not None:
            if not isinstance(float(self.x_min), float):
                raise TypeError("x_min needs to be given a float")
            self.x_min = float(self.x_min)
        if self.y_min is not None:
            if not isinstance(float(self.y_min), float):
                raise TypeError("y_min needs to be given a float")
            self.y_min = float(self.y_min)
        
        if self.x_max is not None:
            if not isinstance(float(self.x_max), float):
                raise TypeError("x_max needs to be given a float")
            self.x_max = float(self.x_max)
        if self.y_max is not None:
            if not isinstance(float(self.y_max), float):
                raise TypeError("y_max needs to be given a float")
            self.y_max = float(self.y_max)
