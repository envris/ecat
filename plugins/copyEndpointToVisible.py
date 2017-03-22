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


from  plugin import reportProvider
from PyQt4 import QtCore,QtGui
from pprint import pprint
import re
from reportDialog import reportDialog
import datetime
from lxml import etree
import sys



class copyDatapointToVisible(reportProvider):
    """Report provider wrapper for the ECOSAR widget."""

    def __init__(self,inputTuple):
        super(copyDatapointToVisible,self).__init__(inputTuple)
        self.reportType = 'Copy data to visible Chemicals ...'


    def displayReport(self):
        """Display report."""
        
        #prompt for save, then close.
        sureMsg = QtGui.QMessageBox()
        sureMsg.setText('Are you sure you want to copy the currently selected data point to all visible chemicals?')
        sureMsg.setStandardButtons(QtGui.QMessageBox.Ok | \
                QtGui.QMessageBox.Cancel)
        sureMsgOut = sureMsg.exec_()

        if sureMsgOut == QtGui.QMessageBox.Ok:
            selectedTi = self.parent.cv.dataPointSelectionModel.currentIndex()
            self.parent.cv.copyDataPointToVisible(selectedTi)
            return
        else:
            return 
