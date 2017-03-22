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


from PyQt4 import QtGui,QtCore,QtWebKit

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


class textDialog(QtGui.QDialog):
    def __init__(self,parent,inputText):
        super(textDialog,self).__init__(parent)
        self.inputText = inputText

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
        self.te = QtGui.QTextEdit()
        #self.wv.setHtml('<html><body><h2>hello!</h2></body></html>')
        #print self.dataLocation
        self.te.setPlainText(self.inputText)
        mFont = QtGui.QFont('Monospace')
        self.te.setFont(mFont)

        vbox.addWidget(self.te)
        vbox.addLayout(hbox)
        self.setLayout(vbox)
        #self.showFullScreen()
        self.setGeometry(50,50,1024,720)
        #self.setGeometry(QtGui.QDesktopWidget().availableGeometry())
        self.show()


    def saveHtml(self):
        saveFileName = QtGui.QFileDialog.getSaveFileName(self,'Save text ...')

        try:
            with open(saveFileName,'w') as saveFile:
                saveFile.write(self.te.doucument().toPlainText())
            print 'Saved as '+saveFileName

        except IOError:
                QtGui.QMessageBox.information(self,'Cascading Viewer',\
                    'Couldn\'t save for some reason.',\
                    QtGui.QMessageBox.Ok)


