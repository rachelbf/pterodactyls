__shortversion__= u'0.1'
__version__= u'0.1.0.dev0' # [N!]N(.N)*[{a|b|rc}N][.postN][.devN]

__all__ = ['TICobject', 'helpers', 'Injection', 'plotter','lightcurve']

from . import helpers, plotter
from .pterodactyls import dactyl
from .injection import Injection
from .lightcurve import LightCurve
from .recovery import IRgrid