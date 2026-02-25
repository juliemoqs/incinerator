#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
from importlib import metadata

#configure the logging
logger = logging.getLogger(__name__)

from .localize import Localize # noqa

__all__ = ["Localize"]

#reads package version from pyproject.toml file
__version__ = metadata.version("incinerator")

#add any global variables
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))