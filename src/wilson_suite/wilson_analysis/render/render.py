from .matplotlib_renderer import MatplotlibRenderer
import numpy as np

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ...wilson_main.spectrum_abstractions import SpecEvalSetup

def render_spectrum(spec_data, spec_eval_setup: 'SpecEvalSetup', 
                    do_diagn, features=None) -> None:
    
    """
    High-level function to render spectrum with specified backend
 -------------
    spec_data is not the final Z values of spectrum (intensities)
        it is amplitudes here, complex values, but they get to be transformed
        using spec_data_operations attribute of RenderingInfo instance

    context = {'spec_data': self.spec, 'system': self.system, 
                'experiment': self.exp, 'diagn': self.diagn, 
                'name': self.name, 
                'spec_eval_setup': self.spec_eval_setup,
                'do_diagn': True}
--------------

    what are the options for rendering? 
    
    by plot type/projection type
    1. 1D plot : slice of a 2D spectrum; 1D IR/Raman spectrum (?)
    2. 2D plot : 2D IR/Raman spectrum; slice of a 3D spectrum : a) contour plot; b) scatter plot
    3. 3D plot : 2D IR/Raman spectrum : surface plot == contout plot (same info; Z axis is intensity or color is intensity)

    by spectrum dimensionality
    - 1D spectrum: 1D IR/Raman spectrum - simply 1d plot
    - 2D spectrum: 2D IR/Raman spectrum - 2D contour/scatter plot or 3D surface plot
    - 3D spectrum: 3D IR/Raman spectrum - 1) 2D slice of a 3D spectrum - 2D contour/scatter plot or 3D surface plot; 2) 3D surface plot with color as intensity
    - nD spectrum: nD IR/Raman spectrum - lower D slices as above
    
    """    
    if not isinstance(spec_data, np.ndarray):
        raise TypeError("spec_data should be a np.ndarray")
    
    if not hasattr(spec_eval_setup, 'is_ready_render'):
        raise TypeError('spec_eval_setup should be a SpecEvalSetup instance')
    
    if not spec_eval_setup.is_ready_render:
        raise ValueError('spec_eval_setup does not have all rendering configs')
    
    if spec_data.size == 0:
        raise ValueError('Empty spec_data array')
    
    if len(spec_data.shape) != 2:
        raise NotImplementedError('only 2D contour plots can be made - input spectrum data is not 2D')

    filename = spec_eval_setup.rnd_info.filename
    backend = spec_eval_setup.rnd_info.backend

    if backend == 'matplotlib':
        renderer_class=MatplotlibRenderer
    else:
        raise NotImplementedError('Only matplotlib backend is currently supported')

    renderer = renderer_class(spec_data=spec_data, 
                              spec_grid=spec_eval_setup.grid,
                              ev_info=spec_eval_setup.ev_info, 
                              rnd_info=spec_eval_setup.rnd_info,
                              do_diagn=do_diagn)
    fig, ax, contour, cbar = renderer.render(filename)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = f'max(intensity) {np.max(np.abs(spec_data)**2):.3e}'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    
    renderer.save_plot(plot_obj=(fig, ax, cbar), filename=filename)

    if features is not None:
        colors = {'a+b,a': 'cyan', 'b,a': 'orange'}
        print('hello')
        for f in features:
            a = f.term_contributions[0].states_parameters[0]['a']
            b = f.term_contributions[0].states_parameters[0]['b']
            res_type = f.term_contributions[0].res_motif.to_str()
            color = colors.get(res_type, 'white')
            A = f.location['A']
            B = f.location['B']
            amp = f.amplitude_coeff
            ax.scatter(A, B, color=color, s=abs(amp)*200,
                    alpha=1.0, edgecolors='white', linewidths=1.5, zorder=5)
            ax.annotate(f"({a},{b})", (A, B),
                        fontsize=7, color='white',
                        xytext=(5, 5), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5))

        # resave with overlaid points
        fig.savefig(filename.replace('.', '_sticks.'), bbox_inches='tight', 
                    dpi=spec_eval_setup.rnd_info.style_config.dpi)
        
    if do_diagn:
        return tuple([fig, ax, contour, cbar]), {'renderer': renderer}
    else:
        return tuple([fig, ax, contour, cbar]), {}
