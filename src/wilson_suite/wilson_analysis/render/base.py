from typing import Any
from typing import TYPE_CHECKING, Protocol
if TYPE_CHECKING:
    from wilson_suite.wilson_main.abstractions import EvaluatedResult, SealedSetup

class Renderer(Protocol):
    """Protocol for a SimulationRun renderer.
    
    A renderer is constructed from (result, sealed) and exposes methods
    for producing visualizations. The protocol only requires the constructor
    signature; view methods are conventions (see MatplotlibRenderer).
    
    Conventional view methods, each returning an Axes-like object for 
    composition (pass `ax` to overlay on an existing figure):
    
        contour(slice_spec=None, ax=None, **opts)
        line(dim, at, ax=None, **opts)
        scatter_features(features=None, label=False, ax=None)
    
    Conventional composites:
    
        contour_with_features(...)
    
    Renderers may implement a subset; users who call a missing method 
    get an AttributeError.
    """
    def __init__(self, result: 'EvaluatedResult', sealed: 'SealedSetup') -> None: ...
    def contour(self, *, slice_spec=None, ax: Any = None, **opts) -> Any: ...
    def line(self, *, dim, at, ax: Any = None, **opts) -> Any: ...
    def scatter_features(self, *, features=None, label=False, ax: Any = None) -> Any: ...



class RendererTemplate:
    """Copy-and-adapt template for writing a custom renderer.
    
    This isn't a base class to inherit — it's a reference. Copy this file,
    rename the class, fill in the methods for your backend. The conventional
    method shapes are stubbed with NotImplementedError so you can see what
    SimulationRun users might call.
    """
    def __init__(self, result, sealed):
        self._result = result
        self._sealed = sealed
    
    def contour(self, *, slice_spec=None, ax=None, **opts):
        raise NotImplementedError
    
    def line(self, *, dim, at, ax=None, **opts):
        raise NotImplementedError
    
    def scatter_features(self, *, features=None, label=False, ax=None):
        raise NotImplementedError
    
    def contour_with_features(self, **opts):
        raise NotImplementedError