"""
POpSiCLE viewer library

This python library contains modules which helps to post-process 
and plot the spectral information obtained with the main routines of
the library.
"""
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from . import constants as const
from . import input_reader as inp
from . import axes as axes
from . import data as data
from . import coords as coords
from . import aux_funcs as aux_funcs
from . import cross as cross
from . import graph as graph


__all__ = ["constants", "input_reader", "axes",
           "data", "coords", "aux_funcs", "cross",
           "graph"]
