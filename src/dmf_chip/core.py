# coding: utf-8
from __future__ import absolute_import, unicode_literals, print_function
import pint

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = ['__version__', 'ureg']

ureg = pint.UnitRegistry()
