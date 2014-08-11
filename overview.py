#!/usr/bin/env python
"""
Gutierrez, B.  Generates the overview.rst file
"""
import numpy as np

def overview_rst:
	file = open("diagnostics.rst","w")
	file.write('Preliminary Diagnostics\n' + '======================\n\n')
	file.write('.. image:: .. /png_dump/medn_tsnr\n')
	file.write('	:width: 49%\n')
	file.write('.. image:: .. /png_dump/tsnr_ratio\n')
	file.write('	:width: 49%\n')
