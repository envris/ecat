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
from lxml import etree
from makeTags import makeTags
from copy import deepcopy
import printElement


class fieldCollection(QtCore.QAbstractListModel):
    """Model containing a list of fields for display."""

    def __init__(self,fieldList=[],parent=None):
        """Constructor.
        
        Takes casrnList (list of ints) and saves it as self.casrnList."""
        super(fieldCollection,self).__init__(parent)

        #This has structure [['Display string','tag string'], ...]
        self.fieldList = fieldList


    def data(self,index,role=QtCore.Qt.DisplayRole):
        """Returns field at index."""
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            #return display string
            #print [self.fieldList[index.row()][0]]
            return self.fieldList[index.row()][0]
        elif role == QtCore.Qt.UserRole:
            #return tag string
            #print [self.fieldList[index.row()][1]]
            return self.fieldList[index.row()][1]
        #elif role == QtCore.Qt.EditRole:
        #    return self.fieldList[index.row()]
        
        return None

    def setData(self,index,value,role=QtCore.Qt.EditRole):
        """Sets data at index with role to value."""
        print 'setData: index, value, role: ',index,value,role
        value = value.toString()
        if role == QtCore.Qt.DisplayRole:
            self.fieldList[index.row()][0]=value
            return True
        elif role == QtCore.Qt.UserRole:
            self.fieldList[index.row()][1]=value
            return True
        else:
            return False

    def itemData(self,index):
        """Gets the data in all roles at index and puts it into a dict."""
        #print 'getItemData'
        variantDict = {}
        for i in range(QtCore.Qt.UserRole+1):
            variantDict[i]=self.data(index,i)
        #print 'getItemData'
        #pprint(variantDict)
        #print 'getItemData done'
        return variantDict

    def setItemData(self,index,variantDict):
        """sets the data at index for all roles defined in variantDict."""
        #print 'setItemData'
        for i in variantDict.keys():
            self.setData(index,variantDict[i],i)
        #pprint(variantDict)
        #print 'setItemData done'
        return True




    def rowCount(self,parent=QtCore.QModelIndex()):
        """Returns length of self.fieldList."""

        if parent.isValid():
            return 0
        else:
            return len(self.fieldList)

    def appendRow(self,newField,parent=QtCore.QModelIndex()):
        """Tacks a field onto the end of self.casrnList."""
        newField[1] = newField[1].toPyObject()
        #print 'appendRow',[newField]

        self.beginInsertRows(parent,len(self.fieldList),len(self.fieldList))
        self.fieldList.append(newField)
        self.endInsertRows()

    def removeRow(self,row,parent=QtCore.QModelIndex()):
        """Removes the field at a particular position."""

        self.beginRemoveRows(parent,row,row)
        self.fieldList.pop(row)
        self.endRemoveRows()

    def insertRows(self,row,count,parent=QtCore.QModelIndex()):
        """Inserts blank rows into the model."""
        self.beginInsertRows(parent,row,row+count-1)
        for idx in range(count):
            self.fieldList.insert(row,[None,None])
        self.endInsertRows()
        return True



    def removeRows(self,row,count,parent=QtCore.QModelIndex()):
        """Removes count fields at a position row."""
        try:
            self.beginRemoveRows(parent,row,row+count-1)
            self.fieldList = self.fieldList[:row] + self.fieldList[(row+count):]
            self.endRemoveRows()

            return True
        except:
            return False


    def supportedDropActions(self):
        """Required for drag/drop functionality."""
        return QtCore.Qt.CopyAction|QtCore.Qt.MoveAction

    def flags(self,index):
        """Returns flags.
        
        This needed to be reimplemented to get the drag/drop functionality
        to work."""
        if index.isValid():
            return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|\
                    QtCore.Qt.ItemIsEnabled
        
        return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|\
                    QtCore.Qt.ItemIsDropEnabled|QtCore.Qt.ItemIsEnabled
    def returnDataList(self):
        """Returns output list of data points."""
        returnList =[]
        for idx in range(self.rowCount()):
            newIdx = self.index(idx)
            returnList.append([self.data(newIdx,QtCore.Qt.DisplayRole),\
                    self.data(newIdx,QtCore.Qt.UserRole)])
        return returnList

#class moveListModel(QtGui.QStringListModel):
#    def __init__(self):
#        super(moveListModel,self).__init__()
#
#
#    def supportedDropActions(self):
#        return QtCore.Qt.CopyAction|QtCore.Qt.MoveAction
#
#    def flags(self,index):
#        if index.isValid():
#            return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|\
#                    QtCore.Qt.ItemIsEnabled
#        
#        return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|\
#                    QtCore.Qt.ItemIsDropEnabled|QtCore.Qt.ItemIsEnabled
#
#    def returnDataList(self):
#        returnList =[]
#        for idx in range(self.rowCount()):
#            print idx
#            newIdx = self.index(idx)
#            returnList.append([self.data(newIdx,QtCore.Qt.DisplayRole).toPyObject(),\
#                    self.data(newIdx,QtCore.Qt.UserRole).toPyObject()])
#        return returnList

class dataPointSelectionDialog(QtGui.QDialog):
    """Data point selection dialog."""
    def __init__(self,parent,model,chemIdx,\
            dpIdx,qIdx,visible):
        super(dataPointSelectionDialog,self).__init__(parent)
        self.parent = parent
        self.model = model
        self.chemIdx = chemIdx
        self.dpIdx = dpIdx
        self.qIdx = qIdx
        self.chemProxyModel = self.parent.cv.chemProxyModel
        self.visible = visible

        self.initUI()


    def initUI(self):
        """Initialise user interface."""
        doneButton = QtGui.QPushButton('Cancel')
        doneButton.clicked.connect(self.accept)

        displayReportButton = QtGui.QPushButton('Display Report')
        displayReportButton.clicked.connect(self.display)


        addButton = QtGui.QPushButton('Add')
        addButton.clicked.connect(self.addField)

        initialFields = [[QtCore.QString(u'CAS-RN'), \
                QtCore.QString(u'CASNUMBER')], \
                [QtCore.QString(u'Common name'), \
                QtCore.QString(u'CommonName')], \
                [QtCore.QString(u'AICS Name'), \
                QtCore.QString(u'AICSName')], \
                [QtCore.QString(u'Assessment Status'), \
                QtCore.QString(u'Status')]]

        #hbox = QtGui.QHBoxLayout()
        #hbox.addStretch(1)
        #hbox.addWidget(addButton)


        #Make a combo box with the contents of the schema in it.
        self.fieldCombo = QtGui.QComboBox()
        for key in sorted(self.model.structureDict.keys(),key = lambda a: \
                self.model.structureDict[a]['longName'].lower()):
            self.fieldCombo.addItem(self.model.structureDict[key]\
                    ['longName'],key)

        self.fieldMdl = fieldCollection(initialFields)
        #self.fieldMdl = moveListModel()

        #QtGui.QAbstractItemView.AboveItem

                
        self.fieldView = QtGui.QListView()
        self.fieldView.setModel(self.fieldMdl)
        self.fieldView.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.fieldView.setDragEnabled(True)
        self.fieldView.setAcceptDrops(True)
        self.fieldView.setDropIndicatorShown(True)


        self.fieldView.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.fieldView.setDragDropMode( QtGui.QAbstractItemView.InternalMove)
        self.fieldView.setDefaultDropAction( QtCore.Qt.MoveAction)
        self.fieldView.setDragDropOverwriteMode(False)

        self.makeContextMenu()
        
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.fieldCombo)
        hbox.addWidget(addButton)

        hbox2 = QtGui.QHBoxLayout()
        hbox2.addStretch(1)
        hbox2.addWidget(doneButton)
        hbox2.addWidget(displayReportButton)


        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(hbox)
        vbox.addWidget(self.fieldView)
        vbox.addLayout(hbox2)
        #vbox.addStretch(1)
        self.setLayout(vbox)

        self.fieldCombo.activated.connect(self.addField)

    def addField(self):
        """Adds a field to the model."""
        rowCount =self.fieldMdl.rowCount()
        idx = self.fieldCombo.currentIndex()
        text = self.fieldCombo.itemText(idx)
        data = self.fieldCombo.itemData(idx)
        self.fieldMdl.appendRow([text,data])
        
        #self.fieldMdl.insertRows(rowCount,1)
        #newIdx = self.fieldMdl.index(rowCount)
        #self.fieldMdl.setData(newIdx,text,QtCore.Qt.DisplayRole)
        #self.fieldMdl.setData(newIdx,data,QtCore.Qt.UserRole)
        #print self.fieldMdl.data(newIdx,QtCore.Qt.UserRole).toPyObject()
        #print [text,data.toPyObject()]

    def makeContextMenu(self):
        '''Slot to make context menu.
        
        Appears on right click in self.fieldView QListView.'''
        self.lvContextMenu = QtGui.QMenu(self)
        self.lvContextMenu.addAction('Delete',self.onActionDeleteItem)
        self.fieldView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.fieldView.customContextMenuRequested.connect(self.onViewContextMenu)
        

    def onViewContextMenu(self,point):
        '''Slot which makes the context menu for self.fieldView.'''
        idx = self.fieldView.indexAt(point)
        if idx.row() >=0:
            self.lvContextMenu.exec_(self.fieldView.mapToGlobal(point))

    def onActionDeleteItem(self):
        '''Action to delete selected item in self.fieldMdl.'''
        selectionModel = self.fieldView.selectionModel()
        idx = selectionModel.currentIndex()

        self.fieldMdl.removeRow(idx.row())

    def display(self):
        """Make and display the group report."""
        def getPosition(qual):
            param = self.model.structureDict[qual.getparent().tag]['qualifiers'][qual.tag]['longName']
            if param in orderList:
                return orderList.index(param)
            else:
                return len(orderList) + 1

        print 'Display!'
        orderList = ['Trophic Level','Species','Acute or chronic endpoint?',\
            'Endpoint name',  'Test duration', 'Pivotal endpoint?', \
            'Measurement type (Measured, Calculated, Read-across, Analogue)', \
            'Comment', 'Citation']

        #get list of stuff to display in the report
        fieldList = self.fieldMdl.returnDataList()
        fieldList = [[str(x[0]),str(x[1])] for x in fieldList]
        #pprint(fieldList)

        html = '''<html>
<head>
<style>
  *{
    font-family: sans-serif;
    }

  .chem {
    position: relative;
    display: inline-block;
    width: 30%;
    /*float: left;*/
    padding: 12px;
    /*border-color : black;
    border-style : solid;
    border-width : 1px;*/
    }

/*
  .contains-three{
    border-color : black;
    border-style : solid;
    border-width : 1px;
  }
*/

  table td{
        -ms-word-break: break-all;
        /*word-break: break-all;*/

        /* Non standard for webkit */
        word-break: break-word;

        -webkit-hyphens: auto;
        -moz-hyphens: auto;
        -ms-hyphens: auto;
            hyphens: auto;

  }
  .word-wrap{
        -ms-word-break: break-all;
        word-break: break-all;

        /* Non standard for webkit */
        word-break: break-word;

        -webkit-hyphens: auto;
        -moz-hyphens: auto;
        -ms-hyphens: auto;
            hyphens: auto;

  }

/* The thing with the tooltip*/
.hover-container {
        position: relative;
        /*border: 1px black solid;*/
        cursor: pointer;
}

/*the tooltip itself*/
.hover-container div {
    position: absolute;
    font-size: 75%;
    display: none;
    background: #000000;
    color: #FFFFFF;
    top: -10px;
    left:   100%;
    border: 1px black solid;
    border-radius: 6px;
    margin-left: 15px;
    padding: 10px;
    opacity: 1.0;
    z-index: 2000;
    /*cross-browser min-width*/
    width: auto !important;
    width:175px;
    min-width:175px;
}
/*The arrow.*/
.hover-container div:after {
  content: '';
  position: absolute;
  top: 10px;
  left: -8px;
  margin-left: -8px;
  width: 0; height: 0;
  opacity: 1.0;
  z-index: 2001;
  border-top: 8px solid transparent;
  border-left: 8px solid transparent;
  border-right: 8px solid #000000;
  border-bottom: 8px solid transparent;
}
/*Display the tooltip*/
.hover-container:hover div{
    display: block;
}
  

</style>
</head>
<body>
'''
        i = -1
        for ti in self.visible:
            i=i+1
            if i%3 == 0:
                html = html + '\n<div class="contains-three">'

            chemXml = ti.getXmlElement()
            svg = ti.getSvg()
            
            if svg == None:
                svg = ''
            else:
                svg = svg+'</br>\n'

            html = html + '<div class="chem word-wrap">\n' +svg
            #html = html+'<table class="word-wrap">\n'


            #printList = []
            for fl in fieldList:
                #1st element in fl is name to print on screen,
                #second element is the name in the XML.
                hits = chemXml.xpath('./'+fl[1]+'[ACTIVE[last()]/Q_QUAL_VALUE1=\'Y\']')
                for index,hit in enumerate(hits):
                    toolTipList = []
                    #get all qualifiers (except ACTIVE)
                    quals = hit.xpath('./*[* and not(self::ACTIVE)]')
                    #print quals
                    for qual in sorted(quals, key=getPosition):
                        longQName = '<b>'+self.model.structureDict[hit.tag]['qualifiers'][qual.tag]['longName']+'</b> '

                        toolTipList.append('<p>'+longQName+printElement.printElement(qual)+'</p>')
                    toolTip = '\n'.join(toolTipList)
                    #pprint(toolTipList)

                    if len(quals)>0:
                        hits[index] = '<span class="hover-container"><u>'+printElement.printElement(hit)+'</u><div>'+toolTip+'</div></span>'
                    else:
                        hits[index] = printElement.printElement(hit)



                #printList.append([fl[0],hits])

                #html = html + '<tr><td><b>'+fl[0]+'</b> '+', '.join(hits)+'</td></tr>\n'
                html = html + '<b>'+fl[0]+'</b> '+', '.join(hits)+'<br/>\n'
            
            #html = html + '\n</table>\n</div>'
            html = html + '\n</div>'

            if i%3 == 2:
                html = html + '\n</div>\n'


        #pprint(printList)
        #print html
        html = html +'\n</body>\n</html>'
        reportDialog(self,html)

                



        return
 



class groupReport(reportProvider):
    """Group report plugin."""

    def __init__(self,inputTuple):
        super(groupReport,self).__init__(inputTuple)
        self.reportType = 'Group report'


    def displayReport(self):
        """Runs dialog to generate group report."""

       
        #get the visible TreeItems
        self.visible = self.getVisible()

        #the TreeItem class has getSvg and getXmlElement
        #methods which return the svg and the etree element.
        #for item in self.visible:
        #    print item.getSvg()
        #    print item.getXmlElement()


        td = dataPointSelectionDialog(self.parent,self.model,self.chemIdx,\
                self.dpIdx,self.qIdx,self.visible)
        dataPointSelectionDialogStatus = td.exec_()
        #if dataPointSelectionDialogStatus:
            #do stuff
