from .core import *
import dmf_chip.core
try:
    from .plot import *
    from .load import *
    import dmf_chip.plot
    import dmf_chip.load
except ImportError:
    import logging

    logging.getLogger(name=__name__).warning('Plotting functions disabled. '
                                             'Please install `matplotlib`.')
