# Copyright (c) 2016 Commonwealth of Australia as represented by
# the Department of the Environment and Energy
#
# This file is free software: you may copy, redistribute and/or modify it
# under the terms of the GPL v2 as published by the Free Software Foundation.
#
# This file is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GPL v2 for more details.
#
# You should have received a copy of the GPL v2 along with this program.
# If not, see <http://www.gnu.org/licenses/gpl-2.0.html>.


import os
import glob

modules = glob.glob(os.path.dirname(__file__)+'/*.py')
#from plugins import * will import all modules except those starting
#with a _ .
__all__ = [os.path.basename(aFileNameThatWillNeverBeUsed)[:-3] for aFileNameThatWillNeverBeUsed in modules if not os.path.basename(aFileNameThatWillNeverBeUsed).startswith('_')]


#I initially thought that plugin.py needed to be at the end. 
#It works without doing this.
#__all__ = [os.path.basename(f)[:-3] for f in modules if (not os.path.basename(f).startswith('_')) and os.path.basename(f)[:-3] != 'plugin']

#append plugin.py
#__all__.append('plugin')

print __all__
