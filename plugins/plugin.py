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
#
# This file incorporates work belonging to the public domain. 
# See <http://martyalchin.com/2008/jan/10/simple-plugin-framework/> for details.

from PyQt4 import QtGui

class pluginMount(type):
    def __init__(cls,name,bases,attrs):

        if not hasattr(cls,'plugins'):
            #set up plugin list
            cls.plugins = []
        else:
            #add class to plugin list.
            cls.plugins.append(cls)

class reportProvider(object):
    __metaclass__ = pluginMount
    #def __init__(self,request,*args,**kwargs):
    def __init__(self,inputTuple=None):
        # The input tuple needs to have the form
        # (model,chemIdx,dpIdx,qIdx,parent)
        # where the model is a TreeModel, and chemIdx,
        # dpIdx, qIdx are the currently selected indices
        # (in the TreeModel basis, not the proxyModel
        # basis. Parent is the parent widget (assessmentTool)

        self.model = inputTuple[0]
        self.chemIdx = inputTuple[1]
        self.dpIdx = inputTuple[2]
        self.qIdx = inputTuple[3]
        self.parent = inputTuple[4]

        self.reportType = ''

    def returnReportType(self):
        return self.reportType

    def makeReport(self):
        pass

    def returnReport(self):
        return None

    def displayReport(self):
        return None

    def getVisible(self):
        '''Returns a list of the visible treeItems.'''

        sfModel = self.parent.cv.chemProxyModel
        rootIdx = sfModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0)))
        rowCount = sfModel.rowCount(rootIdx)
        visible = []
        for idx in range(rowCount):
            visible.append(sfModel.mapToSource(rootIdx.child(idx,0)).internalPointer())

        return visible

class reportDialog(QtGui.QDialog):
    def __init__(self,parent,inputHtml):
        super(reportDialog,self).__init__(parent)
        self.inputHtml = inputHtml

        self.initUI()
    def initUI(self):
        self.setWindowFlags(QtCore.Qt.Window or QtCore.Qt.WindowMaximizeButtonHint)
        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(10)
        closeBtn = QtGui.QPushButton('Close')
        closeBtn.clicked.connect(self.done)
        saveBtn=QtGui.QPushButton('Save ...')
        saveBtn.clicked.connect(self.saveHtml)
        hbox.addWidget(closeBtn)
        hbox.addWidget(saveBtn)

        vbox = QtGui.QVBoxLayout()
        self.wv = QtWebKit.QWebView()
        #self.wv.setHtml('<html><body><h2>hello!</h2></body></html>')
        #print self.dataLocation
        self.wv.setHtml(self.inputHtml)
        vbox.addWidget(self.wv)
        vbox.addLayout(hbox)
        self.setLayout(vbox)
        #self.showFullScreen()
        self.setGeometry(50,50,1024,720)
        #self.setGeometry(QtGui.QDesktopWidget().availableGeometry())
        self.show()


    def saveHtml(self):
        saveFileName = QtGui.QFileDialog.getSaveFileName(self,'Save html ...')

        try:
            with open(saveFileName,'w') as saveFile:
                saveFile.write(self.inputHtml)
            print 'Saved as '+saveFileName

        except IOError:
                QtGui.QMessageBox.information(self,'Cascading Viewer',\
                    'Couldn\'t save for some reason.',\
                    QtGui.QMessageBox.Ok)

