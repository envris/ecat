## Copyright (c) 2016 Commonwealth of Australia as represented by
## the Department of the Environment and Energy
##
## This file is free software: you may copy, redistribute and/or modify it
## under the terms of the GPL v2 as published by the Free Software Foundation.
##
## This file is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
## FOR A PARTICULAR PURPOSE. See the GPL v2 for more details.
##
## You should have received a copy of the GPL v2 along with this program.
## If not, see <http://www.gnu.org/licenses/gpl-2.0.html>.


from time import sleep
import os.path
import shutil
from PyQt4 import QtCore,QtGui
import copy
from lxml import etree
import time
import weakref

class autoSaveWorker(QtCore.QObject):
    """Implements autosave functionality."""
    #error = QtCore.pyqtSignal()
    finished = QtCore.pyqtSignal()
    saveCurrentXml = QtCore.pyqtSignal()

    def __init__(self,cv,autoSaveDir):
        super(autoSaveWorker,self).__init__()
        print 'Starting autosave ...'


    def save(self):
        """Save function (connected to QTimer)."""
        print 'save '

        #emit signal to save (connected to slot in cv)
        self.saveCurrentXml.emit()


    def run(self):
        """Starts the timer."""
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.save)
        self.timer.start(600000)
        #self.timer.start(5000)
    
    def killAutosaveWorker(self):
        """Kills the autosave worker and deletes the thread."""
        self.timer.stop()
        currentThread = QtCore.QThread.currentThread()
        currentThread.quit()
        self.deleteLater()
        currentThread.deleteLater()
        print 'Thread killed.'


 
