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

from PyQt4 import QtCore,QtGui
import sys
from lxml import etree
import pybel
import re
import sys
from pprint import pprint

class intBlankValidator(QtGui.QIntValidator):
    """Integer validator which allows blanks."""
    def __init__(self):
        super(intBlankValidator,self).__init__()

    def validate(self,inString,pos):
        """Re-implemented validate function."""
        if str(inString) == '':
            return (QtGui.QValidator.Acceptable,pos)
        return super(intBlankValidator,self).validate(inString,pos)


class booleanEvaluator(object):
    """Contains methods to evaluate filters and apply boolean logic.

    To use, create a filterCollection object, create booleanEvaluator using 
    booleanEvaluator(inputString,filterCollectionObject), 
    
    """

    def __init__(self,inputString,mdl):
        """Constructor.
        
        Takes string containing names of filters combined with boolean logic, 
        tokenizes it, and saves in self.tokenizedString. The model mdl passed to 
        the constructor is the filterCollection model (subclass of 
        QtCore.QAbstractListModel) which is associated with the string to 
        be parsed. This must contain all of the filters referenced in inputString.
        
        Tokens in inputString are either single words or phrases in quotation marks.
        Allowed boolean operators are: and And AND or Or OR Not not NOT."""

        super(booleanEvaluator,self).__init__()
        #print 'booleanEvaluator.__init__: inputString:',inputString
        self.tokenizedString = self.tokenize(inputString)
        #print 'booleanEvaluator.__init__: tokenized String: ',self.tokenizedString
        self.mdl = mdl

    def tokenize(self,expr):
        """Takes string expr and tokenizes it. Doesn't evaluate filters."""

        #firstly, split by whitespace
        tokenized = expr.split()
        #print 'tokenize: ',tokenized

        #case when there is no filter.
        if tokenized == []:
            return [True]
    
    
        #start at end of list and go backwards to minimise problems with indices.
        index = len(tokenized)
        while index > 0:
            index -= 1
            #print index,tokenized
    
            #firstly, do closing parentheses
            while re.match('.*\)$',tokenized[index]):
                #print index,tokenized
                if tokenized[index] == ')':
                    break
                else:
                    temp =tokenized.pop(index)
                    tokenized[index:index] = [temp[:-1],')']
    
            if tokenized[index] == ')':
                continue
    
            #now do opening parentheses
            index2 = index
            while re.match('\(',tokenized[index2]):
                if tokenized[index2] == '(':
                    break
                else:
                    temp = tokenized.pop(index2)
                    tokenized[index2:index2] = ['(',temp[1:]]
                    index2 += 1
    
            index = index2
    
            if tokenized[index] == '(':
                continue
    
    
    
            #join up stuff in quotation marks.
            #find end of quoted string
            endQuotationMarkMatch = re.match('.*"$',tokenized[index])
            
            if endQuotationMarkMatch:
                
                #check to see whether the string has a quotation mark at the beginning
                beginQuotationMarkMatch = re.match('"',tokenized[index])
                if beginQuotationMarkMatch:
                    #just strip off the quotation marks - they are superfluous.
                    tokenized[index] = tokenized[index].strip('"')
                    continue
                else:
                    if index == 0:
                        # there is no opening quotation mark
                        return None
                    else:
                        #find the index of the token containing the opening quotation mark
                        index2 =index
                        while index2 > 0:
                            index2 -= 1
                            beginQuotationMarkMatch2 = re.match('"',tokenized[index2])
                            if beginQuotationMarkMatch2:
                                #join up the relevant elements of list.
                                newToken = ' '.join(tokenized[index2:index+1])
                                tokenized[index2:index+1]=[newToken]
                                tokenized[index2] = tokenized[index2].strip('"')
                                #reset index
                                index = index2
                                break
                            else:
                                if index2 == 0:
                                    #there is no opening quotation mark
                                    return None
        #print 'end of Tokenize: ', tokenized    
        return tokenized
    def checkTokenizedString(self,parent):
        excludeSet = set(['and','AND','And','or','OR','Or',\
                'not','NOT','Not',')','(',True,False]) 
        for index in range(len(self.tokenizedString)):
            if self.tokenizedString[index] not in excludeSet:
                if self.tokenizedString[index] not in self.mdl.filterDict.keys():
                    #Pop up a message box here.
                    reply = QtGui.QMessageBox.question(parent,'Message',\
                        'Filter '+self.tokenizedString[index]+\
                        ' does not exist.', QtGui.QMessageBox.Ok)

                    return False
        return True

    
    def evaluateFilters(self,xml):
        """Finds filter names in self.tokenizeString, and replaces them with the result of the filter applied to xml.
        
        Returns the modified tokenized string."""

        tokenized = self.tokenizedString[:]
        #print 'tokenized at beginning of evaluateFilters:',tokenized
        #print 'self.tokenizedString at beginning of evaluateFilters:',self.tokenizedString
        #print etree.tostring(xml,pretty_print=True)
        
        excludeSet = set(['and','AND','And','or','OR','Or','not','NOT','Not',')','(',True,False])
        #for index in range(len(tokenized)):
        #    if tokenized[index] not in excludeSet:
        #        if tokenized[index] in self.mdl.filterDict.keys():
        #            tokenized[index] = self.mdl.filterDict[tokenized[index]].applyFilter(xml)
        #        else:
        #            #print [tokenized[index]]
        #            #print str(tokenized[index])+' not in self.mdl.filterDict.keys().'
        #            #print self.mdl.filterDict
        #            return None
        ##add a ( at the beginning and a ) at the end.
        for index in range(len(tokenized)): 
            if tokenized[index] not in excludeSet:
                tokenized[index] = self.mdl.filterDict[tokenized[index]].applyFilter(xml)
        tokenized.append(')')
        tokenized = ['('] + tokenized
   
        #print 'tokenized at the end of evaluateFilters:',tokenized
        #print 'self.tokenizedString at the end of evaluateFilters:',self.tokenizedString
        
        return tokenized
    
    def recursiveEvalBool(self,evaluatedTokens):
        """Applies boolean logic to tokenized string.
        
        Applied after self.evaluateFilters."""

        ip = evaluatedTokens
        #print 'beginning of recursiveEvalBool: evaluateTokens',ip
        #print 'beginning of recursiveEvalBool: self.tokenizedString',self.tokenizedString
    
        index = len(evaluatedTokens)
        while index > 0:
            index -= 1
    
            if evaluatedTokens[index] == ')':
                #find the opening bracket
                p = 1
                index2 = index
                while index2 > 0:
                    index2 -= 1
                    if evaluatedTokens[index2] == ')':
                        p += 1
                    if evaluatedTokens[index2] == '(':
                        p -= 1
                    if p == 0:
                        break
    
                if p == 0:
                    #opening bracket is at index2 so long as p == 0.
                    newTokenList = evaluatedTokens[index2+1:index]
                    if index+1 < len(evaluatedTokens)-1:
                        evaluatedTokens = evaluatedTokens[:index2]+[self.recursiveEvalBool(newTokenList)] + evaluatedTokens[index+1:]
                    else:
                        evaluatedTokens = evaluatedTokens[:index2]+[self.recursiveEvalBool(newTokenList)]
    
                    index = index2
                else:
                    print 'Couldn\'t find opening bracket.'
                    print ip
                    exit()
        #Have now evaluated everything in brackets. Apply boolean logic.
        #print 
        #print '----'
        
        #Find all nots. Evaluate.
        index = len(evaluatedTokens)
        while index > 0:
            index -= 1
            if evaluatedTokens[index] == 'NOT' or evaluatedTokens[index] == 'Not' or evaluatedTokens[index] == 'not':
                if index == len(evaluatedTokens)-1:
                    #There's nothing to negate.
                    print 'Syntax error: nothing to negate.'
                    print index,evaluatedTokens
                    exit()
    
                if type(evaluatedTokens[index+1]) != bool:
                    print 'Syntax error: can\'t negate non-boolean.'
                    print index,evaluatedTokens
                    exit()
    
                result = not evaluatedTokens[index+1]
                #print evaluatedTokens
                #print 'NOT result:',result
                #print index,result,[evaluatedTokens[index],evaluatedTokens[index+1]]
                evaluatedTokens.pop(index)
                evaluatedTokens[index] = not evaluatedTokens[index]
    
    
        #find all ands. Evaluate
        index = len(evaluatedTokens)
        while index > 0:
            index -= 1
            if evaluatedTokens[index] == 'AND' or evaluatedTokens[index] == 'And' or evaluatedTokens[index] == 'and':
                if index == len(evaluatedTokens)-1 or index == 0:
                    #And needs two operands
                    print 'Syntax error: AND needs two operands.'
                    print index,evaluatedTokens
                    exit()
    
                if type(evaluatedTokens[index+1]) != bool or type(evaluatedTokens[index-1]) != bool:
                    print 'Syntax error: can\'t apply AND to  non-booleans.'
                    print index,evaluatedTokens
                    exit()
                result = evaluatedTokens[index-1] and evaluatedTokens[index+1]
                #print evaluatedTokens
                #print 'AND result:',result
                #print index,result,[evaluatedTokens[index-1],evaluatedTokens[index],evaluatedTokens[index+1]]
                evaluatedTokens.pop(index-1)
                evaluatedTokens.pop(index-1)
                evaluatedTokens.pop(index-1)
                evaluatedTokens[index-1:index-1] = [result]
    
        #find all ors. Evaluate
        index = len(evaluatedTokens)
        while index > 0:
            index -= 1
            if evaluatedTokens[index] == 'OR' or evaluatedTokens[index] == 'Or' or evaluatedTokens[index] == 'or':
                if index == len(evaluatedTokens)-1 or index == 0:
                    #And needs two operands
                    print 'Syntax error: OR needs two operands.'
                    print index,evaluatedTokens
                    exit()
    
                if type(evaluatedTokens[index+1]) != bool or type(evaluatedTokens[index-1]) != bool:
                    print 'Syntax error: can\'t apply OR to  non-booleans.'
                    print index,evaluatedTokens
                    exit()
                result = evaluatedTokens[index-1] or evaluatedTokens[index+1]
                #print evaluatedTokens
                #print 'OR result:',result
                #print index,result,[evaluatedTokens[index-1],evaluatedTokens[index],evaluatedTokens[index+1]]
                evaluatedTokens.pop(index-1)
                evaluatedTokens.pop(index-1)
                evaluatedTokens.pop(index-1)
                evaluatedTokens[index-1:index-1] = [result]
    
    
    
        if evaluatedTokens == [False]:
            #print self.tokenizedString
            return False
        elif evaluatedTokens == [True]:
            #print self.tokenizedString
            return True
        else:
            print 'input:',ip
            print 'output:',evaluatedTokens
            exit()

    def applyToChemical(self,idx):
        """Convenience method wrapping self.evaluateFilters and self.recursiveEvalBool."""

        xml = idx.internalPointer().element
        #print 'beginning applyToXml'
        evaluatedTokens = self.evaluateFilters(xml)
        #print 'evaluatedTokens: ',evaluatedTokens
        evaluated = self.recursiveEvalBool(evaluatedTokens)
        #print self.tokenizedString
        #print 'finished applyToXml'
        return evaluated

    def parseString(self,inString):
        """Takes inString, tokenizes it and saves as self.tokenizedString.
        
        When the string to be parsed has changed, this method can be used
        to reset it, instead of making a new booleanEvaluator object."""

        self.tokenizedString = self.tokenize(inString)
        #print 'parseString: ',self.tokenizedString


class filter(object):
    """Base filter class.
    
    Serves as a template for the other filters."""
    def __init__(self,name):
        """Constructor.
        
        Sets self.name to name."""

        self.name = str(name)

    def getName(self):
        """Return self.name."""
        return self.name

    def applyFilter(self,xml):
        """Apply filter to chemical/datapoint with xml element xml.
        
        This needs to be implemented in subclasses."""
        return None

    
    def getToolTip(self):
        """Returns tooltip.
        
        This is used by filterCollection to pass a toolTip to 
        a listView used to edit the filterCollection.
        
        This should be implemented in subclasses."""

        return self.name

class DPFilter(filter):
    def __init__(self,name):
        super(DPFilter,self).__init__(name)
        self.be = None


    def applyQualifierFilters(self,xml):
        """Apply the qualifier filters to datapoint XML element xml.
        
        The datapoints are dealt with in the subclasses."""
        #print 'DPFilter.applyQualifierFilters'
        #print 'self.be',self.be

        tempMdl = filterCollection()
        tempBe = booleanEvaluator('',tempMdl)
        #print 'self.be == tempMdl',self.be == tempMdl

        if self.be != None and tempBe != self.be:
            #print 'applyQualifierFilters: self.be.tokenizedString: ',self.be.tokenizedString

            evaluatedTokens = self.be.evaluateFilters(xml)
            evaluated = self.be.recursiveEvalBool(evaluatedTokens)

            return evaluated
        else:
            return True

class casrnFilter(DPFilter):
    """Filter which checks whether the casrn for the chemical at idx is in self.casrnList."""

    def __init__(self,name,casrnList):
        """Constructor.
        
        Takes casrnList (list of ints) and saves them as self.casrnList."""
        super(casrnFilter,self).__init__(name)
        self.casrnList = casrnList

    def getToolTip(self):
        """Returns list of casrns, separated by commas."""
        return ', '.join([str(x) for x in self.casrnList])

    def applyFilter(self,xml):
        """Checks to see whether any CASNUMBER xml element has a value in self.casrnList."""
        #get CAS-RN elements
        casrnElements = xml.xpath('./CASNUMBER/DP_INTEGER1')
        #print 'casrns',[x.text for x in casrnElements]
        #print 'self.casrnList',self.casrnList

        for dpint in casrnElements:
            if dpint.text in [str(x) for x in self.casrnList]:
                #note that the dpint element is the DP_INTEGER1 element.
                #Need to pass the parent CASNUMBER element to
                #applyQualifierFilters.
                if self.applyQualifierFilters(dpint.getparent()):
                    return True 

        return False



        #self.smilesXPath = etree.XPath('./SMILES/DP_CLOB1')


class smartsFilter(DPFilter):
    """Filter which checks whether any SMILES associated with the chemical matches a SMARTS string."""

    def __init__(self,name,smarts):
        """Constructor.
        
        Takes a SMARTS string, compiles it and saves it in self.smarts and self.pySmarts."""
        
        super(smartsFilter,self).__init__(name)
        self.smarts = smarts
        self.pySmarts = pybel.Smarts(self.smarts)

    def getToolTip(self):
        """Returns the SMARTS for use as a tooltip."""
        return self.smarts

    def applyFilter(self,xml):
        """Checks to see whether any SMILES associated with chemical at idx matches self.pySmarts."""
        
        smilesElements = xml.xpath('./SMILES/DP_CLOB1')

        for smilesElement in smilesElements:
            try:
                mol = pybel.readstring('smi',smilesElement.text)
            except:
                return False
            if self.pySmarts.findall(mol) != []:
                if self.applyQualifierFilters(smilesElement.getparent()):
                    return True 

        return False


    def getSmarts(self):
        """Returns self.smarts."""
        return self.smarts

class filterCollection(QtCore.QAbstractListModel):
    """Model containing list of filters.
    
    Set up for use with a QListView."""
    def __init__(self,filterList=None,parent=None):
        """Constructor."""
        super(filterCollection,self).__init__(parent)

        #GOTCHA: Note that you can't have filterlist=[] as the 
        #default value. Then the filterlist in every instance
        #will be the same, and if you append to the list in
        #one instance, the other instances see the updated
        #list. The following gets around this problem.

        if filterList == None:
            self.filters = []
        else:
            self.filters = filterList

        self.filterDict = {x.name : x for x in self.filters}


    def data(self,index,role=QtCore.Qt.DisplayRole):
        """Returns data at index for various roles."""
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            #return name of filter object
            return self.filters[index.row()].getName()
        
        if role == QtCore.Qt.UserRole:
            #return filter object
            return self.filters[index.row()]

        if role == QtCore.Qt.ToolTipRole:
            #return name of filter object ()
            return self.filters[index.row()].getToolTip()

        #return None if we get this far.
        return None

    def rowCount(self,parent):
        """Return number of filters."""
        if parent.isValid():
            return 0
        else:
            return len(self.filters)

    def appendRow(self,newFilter,parent=QtCore.QModelIndex()):
        """Add a filter to the end of the filter list."""

        self.beginInsertRows(parent,len(self.filters),len(self.filters))
        self.filters.append(newFilter)
        self.filterDict[self.filters[-1].name] = self.filters[-1]
        self.endInsertRows()
    
    def insertRow(self,row,newFilter,parent=QtCore.QModelIndex()):
        """Insert a filter at a given position in the filter list."""

        self.beginInsertRows(parent,row,row)
        self.filters.insert(row,newFilter)
        self.filterDict[self.filters[row].name] = self.filters[row]
        self.endInsertRows()

    def removeRow(self,row,parent=QtCore.QModelIndex()):
        """Remove the filter at a given position in the filter list."""
        name = self.filters[row].name
        self.beginRemoveRows(parent,row,row)
        self.filters.pop(row)
        del self.filterDict[name]
        self.endRemoveRows()

    def index(self,row,column,parent=QtCore.QModelIndex()):
        """Return the index of a filter in a given position."""
        if not self.hasIndex(row, column, parent):
            return QtCore.QModelIndex()

        return self.createIndex(row,0,self.filters[row])

class casrnCollection(QtCore.QAbstractListModel):
    """Model containing a list of CAS-RNs.
    
    This is used to make CAS-RN filter objects."""

    def __init__(self,casrnList=[],parent=None):
        """Constructor.
        
        Takes casrnList (list of ints) and saves it as self.casrnList."""
        super(casrnCollection,self).__init__(parent)

        self.casrnList = casrnList

    def data(self,index,role=QtCore.Qt.DisplayRole):
        """Returns CAS-RN at index for display purposes."""
        if not index.isValid():
            return None

        if role == QtCore.Qt.DisplayRole:
            #return cas-rn
            return str(self.casrnList[index.row()])
        
        return None

    def rowCount(self,parent):
        """Returns length of self.casrnList."""

        if parent.isValid():
            return 0
        else:
            return len(self.casrnList)

    def appendRow(self,newCasrn,parent=QtCore.QModelIndex()):
        """Tacks a CAS-RN onto the end of self.casrnList."""

        self.beginInsertRows(parent,len(self.casrnList),len(self.casrnList))
        self.casrnList.append(newCasrn)
        self.endInsertRows()

    def removeRow(self,row,parent=QtCore.QModelIndex()):
        """Removes the CAS-RN at a particular position."""

        self.beginRemoveRows(parent,row,row)
        self.casrnList.pop(row)
        self.endRemoveRows()



class filterWidget(QtGui.QWidget):
    """Base filterWidget class. Subclassed by DPFilterWidget and QFilterWidget"""
    def __init__(self,structureDict,parent=None):
        super(filterWidget,self).__init__(parent)

        self.structureDict = structureDict
        sys.stdout.write('filterWidget.__init__()\n')
        self.mdl = filterCollection()
        self.be = booleanEvaluator('',self.mdl)

        self.initData()
        self.initUI()

    def initData(self):
        """Initialises data. """
        pass
        

        
    def initUI(self):
        """Initialises user interface."""

        sys.stdout.write('filterWidget.initUI()\n')
        self.view = QtGui.QListView()
        self.view.setModel(self.mdl)

        self.makeContextMenu()

        fhbox = QtGui.QHBoxLayout()
        self.boolLineEdit = QtGui.QLineEdit()
        self.boolLineEdit.setText('')
        self.boolLineEdit.returnPressed.connect(self.setFilter)


        self.setFilterBtn = QtGui.QPushButton('Run Filter')
        self.setFilterBtn.clicked.connect(self.setFilter)
        fhbox.addWidget(self.boolLineEdit)
        fhbox.addWidget(self.setFilterBtn)


        self.layout = QtGui.QVBoxLayout()
        self.layout.addLayout(fhbox)
        self.layout.addWidget(self.view)

        self.setLayout(self.layout)
        self.setWindowTitle('Testing filter view')


    def setFilter(self):
        """Sets the filter using contents of self.boolLineEdit."""
        try:
            self.be.parseString(str(self.boolLineEdit.text()))
        except:
            reply = QtGui.QMessageBox.information(self,'Message',
                    'Something is wrong with your boolean logic!',QtGui.QMessageBox.Ok)
        #print 'Test SMARTS search:', self.be.applyToXml(self.xml)



    def makeContextMenu(self):
        """Slot to make context menu.
        
        Appears on right click in self.view QListView."""
        self.lvContextMenu = QtGui.QMenu(self)
        self.lvContextMenu.addAction('Delete selected filter',self.onActionDeleteFilter)
        self.view.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.view.customContextMenuRequested.connect(self.onFilterContextMenu)

    def onFilterContextMenu(self,point):
        """Slot which makes context menu for the self.listView."""

        idx = self.view.indexAt(point)
        if idx.row() >=0:
            self.lvContextMenu.exec_(self.view.mapToGlobal(point))

    def onActionDeleteFilter(self):
        """Action to delete selected filter."""
        
        selectionModel = self.view.selectionModel()
        idx = selectionModel.currentIndex()

        self.mdl.removeRow(idx.row())

    def getBooleanEvaluator(self):
        """Return self.be.
        
        Note that this contains a pointer to self.mdl, so a getFilterCollection
        method is not required."""
        return self.be


class DPFilterWidget(filterWidget):
    """This widget displays the available filters and allows the user to add new filters, edit existing ones, and combine them using boolean logic."""
    def __init__(self,structureDict,parent=None):
        """Constructor.
        
        Initialise data and user interface."""
        super(DPFilterWidget,self).__init__(structureDict,parent)


    def initUI(self):
        """Initialises user interface."""
        sys.stdout.write('DPFilterWidget.initUI()\n')
        super(DPFilterWidget,self).initUI()

        hbox = QtGui.QHBoxLayout()
        addFilterBtn = QtGui.QPushButton('Add Filter')

        addFilterBtn.clicked.connect(self.runFilterDialog)

        hbox.addWidget(addFilterBtn)
        self.layout.addLayout(hbox)

        #self.view.doubleClicked.connect(self.determineEditType)
        self.view.doubleClicked.connect(self.edit)

    def edit(self,idx):
        """Edit filter at idx."""
        underlyingObject=idx.internalPointer()
        #print 'DPFilterWidget.edit(self,idx): underlyingObject',underlyingObject
        self.runFilterDialog(underlyingObject,idx)

    def applyFilter(self,idx):
        """Run the boolean evaluator on chemical at index idx.
        
        This index is from the underlying TreeModel used in the cascading viewer.
        The method is used to work out whether the chemical should be visible or
        not."""

        return self.be.applyToChemical(self,idx)

    def runFilterDialog(self,inputFilter=None,idx=None):
        """Run the filter dialog box (instance of DPFilterDialog)."""
        #Fix what QPushButton spits out.
        #print 'DPFilterWidget.runFilterDialog: ',inputFilter
        if inputFilter == False:
            inputFilter = None
        #print 'DPFilterWidget.runFilterDialog: modified',inputFilter
        if inputFilter == None:
            fd = DPFilterDialog(self.structureDict,self)
        else:
            fd = DPFilterDialog(self.structureDict,self,inputFilter)

        fdStatus = fd.exec_()
        if fdStatus:
            #get the new filter
            newFilter = fd.constructFilter()
            #print 'runFilterDialog: newFilter:',newFilter
            #print 'runFilterDialog: currentFilterWidget:',fd.currentFilterWidget
            #print 'type(fd.currentFilterWidget)',type(fd.currentFilterWidget)
            #print 'type(newFilter)',type(newFilter)


            #bail if we couldn't make the new filter for some reason.
            if not issubclass(type(newFilter),filter):
                    msg = QtGui.QMessageBox()
                    msg.setText('Couldn\'t make filter!')
                    msg.exec_()
                    return

            name = newFilter.getName()
            #print 'newFilter.getName: ',[name]
            #print 'current names: ', [str(x.name) for x in self.mdl.filters]

            if inputFilter == None:
                if name not in [str(x.name) for x in self.mdl.filters]:
                    self.mdl.appendRow(newFilter)
                else:
                    msg = QtGui.QMessageBox()
                    msg.setText('Already have a filter with name \"'+name+'\".')
                    msg.exec_()
            else:
                print 'popping old filter'
                #pop the old filter and insert the new one.
                #pop the row in question.
                row = idx.row()
                self.mdl.beginRemoveRows(QtCore.QModelIndex(),row,row)
                tempFilter = self.mdl.filters.pop(row)
                self.mdl.endRemoveRows()

                self.mdl.insertRow(row,newFilter)



#class CLOBQualifierFilter(filter):
#    def __init__(self,name,qualifierName,regexp=None,mostRecent=None,):
#        """Constructor.
#        
#        Note, mostRecent is used to restrict the search to the most recent
#        qualifier of the given type."""
#        
#
#        super(CLOBQualifierFilter,self).__init__(name)
#        self.qualifierName = qualifierName
#        #mostRecent must be either True or False
#        self.mostRecent = mostRecent
#
#        cs = QtCore.Qt.CaseInsensitive
#        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
#        self.regexp=QtCore.QRegExp(regexp,cs,syntax)
#
#
#
#    def getToolTip(self):
#        
#        return str(self.value)
#
#
#    def applyFilter(self,idx):
#        """Checks to see whether any qualifier of type self.qualifierName has a value that matches self.value."""
#        xml = idx.internalPointer().element
#        
#        #get values from xml 
#        values = set([x.text for x in xml.xpath('./'+self.qualifierName+'/Q_CLOB1')])
#        print 'values',values
#        print 'self.value',self.value
#
#        if len(values) > 0:
#            if mostRecent == True:
#                #check most recent value only.
#                if self.regexp.indexIn(values[-1]) >= 0:
#                    return True
#                else:
#                    return False
#            else:
#                #check all values for a match.
#                for v in values:
#                    if self.regexp.indexIn(v) >= 0:
#                        return True
#                return False
#
#        return False
#
#
#
#class Q_VARCHAR2QualifierFilter(filter):
#    def __init__(self,name,qualifierName,regexp=None,mostRecent=None,):
#        """Constructor.
#        
#        Note, mostRecent is used to restrict the search to the most recent
#        qualifier of the given type."""
#        
#
#        super(Q_VARCHAR2QualifierFilter,self).__init__(name)
#        self.qualifierName = qualifierName
#        #mostRecent must be either True or False
#        self.mostRecent = mostRecent
#
#        cs = QtCore.Qt.CaseInsensitive
#        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
#        self.regexp=QtCore.QRegExp(regexp,cs,syntax)
#
#
#
#    def getToolTip(self):
#        
#        return str(self.value)
#
#
#    def applyFilter(self,idx):
#        """Checks to see whether any qualifier of type self.qualifierName has a value that matches self.value."""
#        xml = idx.internalPointer().element
#        
#        #get values from xml 
#        values = set([x.text for x in xml.xpath('./'+self.qualifierName+'/Q_QUAL_VALUE1')])
#        print 'values',values
#        print 'self.value',self.value
#
#        if len(values) > 0:
#            if mostRecent == True:
#                #check most recent value only.
#                if self.regexp.indexIn(values[-1]) >= 0:
#                    return True
#                else:
#                    return False
#            else:
#                #check all values for a match.
#                for v in values:
#                    if self.regexp.indexIn(v) >= 0:
#                        return True
#                return False
#
#        return False

class QFilterWidget(filterWidget):

    def __init__(self,structureDict,parent=None):
        """Constructor.
        
        Initialise data and user interface."""
        sys.stdout.write('QFilterWidget.__init__()\n')
        super(QFilterWidget,self).__init__(structureDict,parent)

    def initData(self):

        #temporarily hard code a QYNFilter for testing.

        #tempFilter = QYNFilter('Active','ACTIVE',False,'Y')
        #self.mdl.appendRow(tempFilter)
        
        #Add an active filter as default.
        activeFilter = QYNFilter('Active','ACTIVE',True,'Y')
        self.mdl.appendRow(activeFilter)




    def initUI(self):
        """Initialises user interface."""
        super(QFilterWidget,self).initUI()

        addFilterBtn = QtGui.QPushButton('Add filter')
        self.layout.addWidget(addFilterBtn)
        addFilterBtn.clicked.connect(self.runQFilterDialog)

        #set active filter as default.
        self.boolLineEdit.setText('Active')
        self.be.parseString(str(self.boolLineEdit.text()))

    
    def runQFilterDialog(self,inputFilter=None,idx=None):
        """Run the filter dialog box (instance of QFilterDialog)."""
        #Fix what QPushButton spits out.
        #print 'QFilterWidget.runQFilterDialog: ',inputFilter
        if inputFilter == False:
            inputFilter = None
        #print 'DPFilterWidget.runFilterDialog: modified',inputFilter
        if inputFilter == None:
            fd = QFilterDialog(self.structureDict,self.currentDpName,self)
        else:
            fd = QFilterDialog(self.structureDict,self.currentDpName,self,inputFilter)

        fdStatus = fd.exec_()
        if fdStatus:
            #get the new filter
            newFilter = fd.constructFilter()
            #print 'runQFilterDialog: newFilter:',newFilter
            #print 'runFilterDialog: currentFilterWidget:',fd.currentFilterWidget
            #print 'type(fd.currentFilterWidget)',type(fd.currentFilterWidget)
            #print 'type(newFilter)',type(newFilter)


            #bail if we couldn't make the new filter for some reason.
            if not issubclass(type(newFilter),filter):
                    msg = QtGui.QMessageBox()
                    msg.setText('Couldn\'t make filter!')
                    msg.exec_()
                    return

            name = newFilter.getName()
            #print 'newFilter.getName: ',[name]
            #print 'current names: ', [str(x.name) for x in self.mdl.filters]

            if inputFilter == None:
                if name not in [str(x.name) for x in self.mdl.filters]:
                    self.mdl.appendRow(newFilter)
                else:
                    msg = QtGui.QMessageBox()
                    msg.setText('Already have a filter with name \"'+name+'\".')
                    msg.exec_()
            else:
                #print 'popping old filter'
                #pop the old filter and insert the new one.
                #pop the row in question.
                row = idx.row()
                self.mdl.beginRemoveRows(QtCore.QModelIndex(),row,row)
                tempFilter = self.mdl.filters.pop(row)
                self.mdl.endRemoveRows()

                self.mdl.insertRow(row,newFilter)

    def setCurrentDpName(self,currentDpName):
        """Sets current dp name.
        
        Called by DPFilterDialog.selectDataType when filterTypeCombo 
        index is changed."""
        self.currentDpName = currentDpName
        #print 'QFilterWidget.currentDpName: ',self.currentDpName

 


class DPFilterDialog(QtGui.QDialog):
    """Dialog to create a new datapoint filter or edit an existing one.
    
    When self.filterTypeCombo is changed, a new widget of type makeFilter is 
    made in leftVBox. When the dialog is accepted by clicking OK, the makeFilter
    widget's makeFilter method is called by DPFilterWidget.runFilterDialog  to 
    actually make the filter."""

    def __init__(self,structureDict,parent=None,inputFilter=None):

        super(DPFilterDialog,self).__init__(parent)
        self.shortFromLongNameDict = { structureDict[key]['longName'] : \
                key for key in structureDict.keys()}
        self.structureDict = structureDict
        #print 'DPFilterDialog: structureDict: ',self.structureDict
        self.inputFilter = inputFilter
        self.parent = parent
        self.initUI()
        self.initData()

    def initData(self):
        """Initialise data."""
        pass


    def initUI(self):
        """Initialise interface."""

        #print 'DPFilterDialog.initData ',self.inputFilter
        #print 'type(DPFilterDialog.initData) ',type(self.inputFilter)
        
        okBtn = QtGui.QPushButton('OK')
        okBtn.clicked.connect(self.accept)

        cancelBtn = QtGui.QPushButton('Cancel')
        cancelBtn.clicked.connect(self.reject)

        self.filterTypeCombo = QtGui.QComboBox()
        #self.filterTypeCombo.currentIndexChanged.connect(self.selectDataType)
        self.filterTypeCombo.activated.connect(self.selectDataType)
        self.qf = QFilterWidget(self.structureDict,self)
        #print 'qfilterwidget model filterCollection pointer',hex(id(qf.mdl))
        #print 'parent model filterCollection  pointer',hex(id(self.parent.mdl))
        #print 'qfilterwidget model filterCollection filters pointer',hex(id(qf.mdl.filters))
        #print 'parent model filterCollection filters pointer',hex(id(self.parent.mdl.filters))

        self.leftVBox = QtGui.QVBoxLayout()
        self.rightVBox = QtGui.QVBoxLayout()
        self.leftVBox.addWidget(self.filterTypeCombo)
        self.rightVBox.addWidget(self.qf)
        self.hbox = QtGui.QHBoxLayout()
        self.hbox.addLayout(self.leftVBox)
        self.hbox.addLayout(self.rightVBox)

        self.hbox2 = QtGui.QHBoxLayout()
        self.hbox2.addStretch(1)
        self.hbox2.addWidget(okBtn)
        self.hbox2.addWidget(cancelBtn)

        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addStretch(1)
        self.vbox.addLayout(self.hbox)
        self.vbox.addLayout(self.hbox2)

        self.setLayout(self.vbox)
        
        if self.inputFilter == None:
            #New datapoint.
            #print 'self.inputFilter == None'

            #Populate combobox. The data must be a subclass of makeFilter
            self.filterTypeCombo.insertItem(0,'CAS-RN filter',makeCasrnFilter)
            self.filterTypeCombo.insertItem(1,'SMARTS filter',makeSmartsFilter)
            #self.filterTypeCombo.insertItem(1,'B',Q_VARCHAR2QualifierFilter)


            self.filterTypeCombo.insertSeparator(2)

            #test:
            for key in sorted(self.structureDict.keys(),key = lambda a: \
                    self.structureDict[a]['longName'].lower()):
                if self.structureDict[key]['dpType'] == 'DP_CLOB':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPCLOBFilter)
                elif self.structureDict[key]['dpType'] == 'DP_VARCHAR2':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPVARCHAR2Filter)
                elif self.structureDict[key]['dpType'] == 'DP_YN':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPYNFilter)
                elif self.structureDict[key]['dpType'] == 'DP_VARCHAR2CHOICE':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPVARCHAR2CHOICEFilter)
                elif self.structureDict[key]['dpType'] == 'DP_FLOAT':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPFloatFilter)
                elif self.structureDict[key]['dpType'] == 'DP_INTEGER':
                    self.filterTypeCombo.addItem(self.structureDict[key]\
                            ['longName'],makeDPINTEGERFilter)
                else:
                    print self.structureDict[key]['dpType'] +' isn\'t one of the accepted types.' 

            self.selectDataType(0)

        else:
            #print 'self.inputFilter != None'
            #make sure we can't change the type of the datapoint.
            if type(self.inputFilter) == casrnFilter:
                self.filterTypeCombo.insertItem(0,'CAS-RN filter',makeCasrnFilter)
                self.removeMakeFilters()

                #note that the first argument (self) is the parent.
                #print 'type(self.inputFilter):', type(self.inputFilter)
                self.currentFilterWidget = makeCasrnFilter(self,self.inputFilter)
                #self.selectDataType(0)
            elif type(self.inputFilter) == smartsFilter:
                self.filterTypeCombo.insertItem(0,'SMARTS filter',makeSmartsFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeSmartsFilter(self,self.inputFilter)

            elif type(self.inputFilter) == DPCLOBFilter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPCLOBFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPCLOBFilter(self,dpName,\
                        longDpName,self.inputFilter)

            elif type(self.inputFilter) == DPVARCHAR2Filter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPVARCHAR2Filter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPVARCHAR2Filter(self,dpName,\
                        longDpName,self.inputFilter)

            elif type(self.inputFilter) == DPYNFilter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPYNFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPYNFilter(self,dpName,\
                        longDpName,self.inputFilter)
            elif type(self.inputFilter) == DPVARCHAR2CHOICEFilter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPVARCHAR2CHOICEFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPVARCHAR2CHOICEFilter(self,dpName,\
                        longDpName,self.structureDict,self.inputFilter)
            elif type(self.inputFilter) == DPFloatFilter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPFloatFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPFloatFilter(self,dpName,\
                        longDpName,self.inputFilter)
            elif type(self.inputFilter) == DPINTEGERFilter:
                dpName = self.inputFilter.dpName
                longDpName = self.inputFilter.longDpName
                self.filterTypeCombo.insertItem(0,longDpName,makeDPINTEGERFilter)
                self.removeMakeFilters()
                self.currentFilterWidget = makeDPINTEGERFilter(self,dpName,\
                        longDpName,self.inputFilter)

            else:
                print str(type(self.inputFilter))+' isn\'t one of the accepted types.' 
                self.reject()

            self.leftVBox.addWidget(self.currentFilterWidget)

            #Now, populate the QFilterWidget with the correct stuff from 
            #self.inputFilter.
            #print 'setting qfilters...'
            self.qf.mdl = self.inputFilter.be.mdl
            boolText = ' '.join(self.inputFilter.be.tokenizedString)
            self.qf.boolLineEdit.setText(boolText)

            #print boolText
            self.qf.be = booleanEvaluator(boolText,self.qf.mdl)
            self.qf.view.setModel(self.qf.mdl)
        
        #self.filterTypeCombo.currentIndexChanged.connect(self.selectDataType)

    def selectDataType(self,idx):
        """Slot called when index of self.filterTypeCombo is changed."""
        #print idx,self.filterTypeCombo.itemData(idx).toPyObject()
        newClass = self.filterTypeCombo.itemData(idx).toPyObject()
        
        if newClass == makeSmartsFilter or newClass == makeCasrnFilter:
            newWidget = newClass(self,self.inputFilter)
            if newClass == makeSmartsFilter:
                self.qf.setCurrentDpName('SMILES')
            elif newClass == makeCasrnFilter:
                self.qf.setCurrentDpName('CASNUMBER')
        else:
            longName = str(self.filterTypeCombo.itemText(idx))
            shortName = self.shortFromLongNameDict[longName]
            if newClass == makeDPVARCHAR2CHOICEFilter:
                #need structureDict
                newWidget = newClass(self,shortName,longName,self.structureDict,self.inputFilter)
            else:
                #need one of the generic makeFilter widgets, which need the name
                #of the datapoint to be passed to them.
                newWidget = newClass(self,shortName,longName,self.inputFilter)
        
            #pass the short name to the filter widget. It needs this to
            #know what qualifiers to put in the QFilterDialog box when
            #the "Add Filter" button is pressed.
            self.qf.setCurrentDpName(shortName)




        self.removeMakeFilters()
        
        self.leftVBox.invalidate()
        self.rightVBox.invalidate()
        self.hbox.invalidate()
        self.hbox2.invalidate()
        self.vbox.invalidate()

        self.currentFilterWidget = newWidget
        self.leftVBox.addWidget(self.currentFilterWidget)

        self.adjustSize()

    def removeMakeFilters(self):
        """Gets rid of makeFilter objects in leftVBox."""
        #get rid of any makeFilter widgets (should be only one widget anyway)
        for index in range(self.leftVBox.count()-1,-1,-1):
            #print type(self.leftVBox.itemAt(index).widget())
            if issubclass(type(self.leftVBox.itemAt(index).widget()),makeFilter):
                item = self.leftVBox.takeAt(index)
                item.widget().deleteLater()

    def constructFilter(self):
        """Constructs filter."""
        newFilter = self.currentFilterWidget.makeFilter()
        if newFilter != None:
        	#make temp filters for comparison.
        	tempMdl = filterCollection()
        	tempBe = booleanEvaluator('',tempMdl)

        	if self.qf.be  != tempBe:
        	    #There is a qualifier filter to add.
        	    newFilter.be = self.qf.be
        	    #print 'DPFilterWidget.constructFilter: ',newFilter.be

        return newFilter

class QFilterDialog(QtGui.QDialog):
    """QFilter dialog class."""

    def __init__(self,structureDict,currentDpName,parent=None,inputFilter=None):

        super(QFilterDialog,self).__init__(parent)
        self.currentDpName = currentDpName
        #print 'QFilterDialog.currentDpName: ',self.currentDpName
        self.structureDict = structureDict
        self.qualifierDict = self.structureDict[currentDpName]['qualifiers']
        #print 'QFilterDialog.qualifierDict',self.qualifierDict
        self.shortFromLongNameDict = {self.qualifierDict[key]['longName'] : \
               key for key in self.qualifierDict.keys()}
        self.inputFilter = inputFilter
        self.parent = parent
        self.initUI()
        self.initData()

    def initData(self):
        """Initialise data."""
        pass

    def initUI(self):
        """Initialise interface."""

        #print 'QFilterDialog.initData ',self.inputFilter
        #print 'type(QFilterDialog.initData) ',type(self.inputFilter)
        
        okBtn = QtGui.QPushButton('OK')
        okBtn.clicked.connect(self.accept)

        cancelBtn = QtGui.QPushButton('Cancel')
        cancelBtn.clicked.connect(self.reject)

        self.filterTypeCombo = QtGui.QComboBox()
        #self.filterTypeCombo.currentIndexChanged.connect(self.selectDataType)
        self.filterTypeCombo.activated.connect(self.selectDataType)

        self.vBox = QtGui.QVBoxLayout()
        self.vBox.addStretch(1)

        self.innerVBox = QtGui.QVBoxLayout()
        self.innerVBox.addWidget(self.filterTypeCombo)

        self.hbox = QtGui.QHBoxLayout()
        self.hbox.addStretch(1)
        self.hbox.addWidget(okBtn)
        self.hbox.addWidget(cancelBtn)

        self.vBox.addLayout(self.innerVBox)
        self.vBox.addLayout(self.hbox)

        self.setLayout(self.vBox)

        qualifierDict = self.structureDict[self.currentDpName]['qualifiers']

        #print qualifierDict
        

        
        if self.inputFilter == None:
            #New datapoint.
            print 'self.inputFilter == None'

            for key in sorted(qualifierDict.keys(),key = lambda a: \
                    qualifierDict[a]['longName'].lower()):
                if qualifierDict[key]['qType'] == 'Q_YN':
                    self.filterTypeCombo.addItem(qualifierDict[key]\
                            ['longName'],makeQYNFilter)
                elif qualifierDict[key]['qType'] == 'Q_CLOB':
                    self.filterTypeCombo.addItem(qualifierDict[key]\
                            ['longName'],makeQCLOBFilter)
                elif qualifierDict[key]['qType'] == 'Q_VARCHAR2':
                    self.filterTypeCombo.addItem(qualifierDict[key]\
                            ['longName'],makeQVARCHAR2Filter)
                else:
                    print key +' isn\'t one of the accepted types.' 

            self.selectDataType(0)


    def selectDataType(self,idx):
        """Slot called when index of self.filterTypeCombo is changed."""
        #print idx,self.filterTypeCombo.itemData(idx).toPyObject()
        newClass = self.filterTypeCombo.itemData(idx).toPyObject()
        
        #need one of the generic makeFilter widgets, which need the name
        #of the datapoint to be passed to them.
        longName = str(self.filterTypeCombo.itemText(idx))
        shortName = self.shortFromLongNameDict[longName]
        newWidget = newClass(self,shortName,longName,self.inputFilter)
        

        self.removeMakeFilters()
        
        self.hbox.invalidate()
        self.innerVBox.invalidate()
        self.vBox.invalidate()

        self.currentFilterWidget = newWidget
        self.innerVBox.addWidget(self.currentFilterWidget)

        self.adjustSize()

    def removeMakeFilters(self):
        """Gets rid of makeFilter objects in leftVBox."""
        #get rid of any makeFilter widgets (should be only one widget anyway)
        #print 'QFilterDialog.removeMakeFilters'
        for index in range(self.innerVBox.count()-1,-1,-1):
            if issubclass(type(self.innerVBox.itemAt(index).widget()),makeFilter):
                item = self.innerVBox.takeAt(index)
                item.widget().deleteLater()

    def constructFilter(self):
        """Constructs filter."""
        newFilter = self.currentFilterWidget.makeFilter()

        return newFilter

        
        
class makeFilter(QtGui.QWidget):
    """Base widget for making filters (for either qualifiers or datapoints)."""
    def __init__(self,parent,inputFilter=None):
        super(makeFilter,self).__init__(parent)

        self.inputFilter = inputFilter

        self.initUI()
        self.initData()

    def initUI(self):
        """Initialise User interface."""
        pass

    def initData(self):
        """Populate widget with data from inputName and inputFilter.""" 
        #implement in subclass
        return None

    def makeFilter(self):
        """Returns filter."""
        #implement in subclass
        return None

class myListView(QtGui.QListView):
    """QListView that accepts lists of CAS-RNs in a drag and drop event."""
    def __init__(self):
        super(myListView,self).__init__()
        self.setAcceptDrops(True)

    def dragEnterEvent(self,e):
            
        if e.mimeData().hasFormat('text/plain'):
            #print 'hasFormat(\'text/plain\')', e.mimeData().text()
            e.accept()
        else:
            e.ignore() 
    def dragMoveEvent(self,e):
        if e.mimeData().hasFormat('text/plain'):
            e.setDropAction(QtCore.Qt.CopyAction)
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        #strip, then add a newline so that the regex is easy to write...
        eventText = str(e.mimeData().text()).strip()+'\n'
        if re.match('^([ \t]*[0-9]+[ \t]*\n)+$',eventText):
            mdl = self.model()
            casrnList = eventText.split()
            for casrn in casrnList:
                if int(casrn) not in set(mdl.casrnList):
                    mdl.appendRow(int(casrn))

        #else do nothing.
        
        


class makeCasrnFilter(makeFilter):
    """Widget to make a new CAS-RN filter or edit an existing one."""
    def __init__(self,parent,inputFilter=None):
        #print 'makeCasrnFilter.parent, inputFilter: ',parent,inputFilter
        super(makeCasrnFilter,self).__init__(parent,inputFilter)

        #Note that superclass runs initUI, then initData.
        #self.inputFilter is set to the contents of inputFilter

    def initData(self):
        """Initialise data with initList (list of CAS-RNs as ints.)"""
        #print 'makeCasrnFilter.inputFilter: ',self.inputFilter
        if self.inputFilter == None:
            self.mdl = casrnCollection([])
            self.nameLineEdit.setText('')
        else:
            self.mdl = casrnCollection(self.inputFilter.casrnList)
            self.nameLineEdit.setText(self.inputFilter.name)

        self.lv.setModel(self.mdl)




    def initUI(self):
        """Initialise user interface."""

        #self.lv = QtGui.QListView()
        self.lv = myListView()
        #self.lv.setModel(self.mdl)

        self.makeContextMenu()

        self.nameLineEdit = QtGui.QLineEdit()
        #self.nameLineEdit.setText(name)
        self.casrnLineEdit = QtGui.QLineEdit()
        iv = QtGui.QIntValidator()
        self.casrnLineEdit.setValidator(iv)
        self.casrnLineEdit.setDragEnabled(True)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)

        hblo = QtGui.QHBoxLayout()
        lbl = QtGui.QLabel('CAS-RN')
        hblo.addWidget(lbl)
        hblo.addWidget(self.casrnLineEdit)
        addBtn = QtGui.QPushButton('Add to list')
        addBtn.clicked.connect(self.addCASRN)
        self.casrnLineEdit.returnPressed.connect(self.addCASRN)
        hblo.addWidget(addBtn)



        #okBtn = QtGui.QPushButton('OK')
        #cancelBtn = QtGui.QPushButton('Cancel')
        #okBtn.clicked.connect(self.accept)
        #cancelBtn.clicked.connect(self.reject)

        vblo = QtGui.QVBoxLayout()
        hblo2 = QtGui.QHBoxLayout()
        #hblo2.addWidget(okBtn)
        #hblo2.addWidget(cancelBtn)

        vblo.addLayout(flo)
        vblo.addWidget(self.lv)
        vblo.addLayout(hblo)
        vblo.addWidget(addBtn)
        vblo.addLayout(hblo2)


        self.setLayout(vblo)
        self.nameLineEdit.setFocus()

    def makeContextMenu(self):
        """Slot to make context menu.
        
        Appears on right click in self.view QListView."""

        self.lvContextMenu = QtGui.QMenu(self)
        self.lvContextMenu.addAction('Delete selected CAS-RN',self.onActionDeleteFilter)
        self.lv.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.lv.customContextMenuRequested.connect(self.onFilterContextMenu)

    def onFilterContextMenu(self,point):
        """Slot which makes context menu for the listview."""

        idx = self.lv.indexAt(point)
        if idx.row() >=0:
            self.lvContextMenu.exec_(self.lv.mapToGlobal(point))

    def onActionDeleteFilter(self):
        """Deletes selected filter."""

        selectionModel = self.lv.selectionModel()
        idx = selectionModel.currentIndex()

        self.mdl.removeRow(idx.row())


    def addCASRN(self):
        """Add CAS-RN in self.casrnLineEit to casrnCollection self.mdl."""
        if self.casrnLineEdit.text() != '':
            #validator makes sure self.casrnLineEdit contains
            #an int or ''
            if int(self.casrnLineEdit.text()) not in set(self.mdl.casrnList):
                self.mdl.appendRow(int(self.casrnLineEdit.text()))
                self.casrnLineEdit.setText('')

    def makeFilter(self):
        """Makes the filter."""
        try:
            returnFilter = casrnFilter(self.nameLineEdit.text(),self.mdl.casrnList)
        except:
            returnFilter = None
        return returnFilter

class makeSmartsFilter(makeFilter):
    """Makes a SMARTS filter."""
    def __init__(self,parent,inputFilter=None):
        """Constructor."""
        #print 'makeSmartsFilter.parent, inputFilter: ',parent,inputFilter
        super(makeSmartsFilter,self).__init__(parent,inputFilter)
        #print 'makeSmartsFilter.inputFilter',self.inputFilter

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initSmarts = self.inputFilter.getSmarts()
        else:
            name = ''
            initSmarts = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.smartsLineEdit =  QtGui.QLineEdit()
        self.smartsLineEdit.setText(initSmarts)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('SMARTS',self.smartsLineEdit)


        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile smarts to see it's valid.
            smarts = str(self.smartsLineEdit.text())
            pySmarts = pybel.Smarts(smarts)
            returnFilter = smartsFilter(name,smarts)

        except:
            returnFilter = None
            QtGui.QMessageBox.information(self,'SMARTS filter',\
                        'SMARTS invalid.',\
                        QtGui.QMessageBox.Ok)

        return returnFilter



class makeDPCLOBFilter(makeFilter):
    """Makes a SMARTS filter."""
    def __init__(self,parent,dpName,longDpName,inputFilter=None):
        """Constructor."""
        self.dpName = dpName
        self.longDpName = longDpName
        super(makeDPCLOBFilter,self).__init__(parent,inputFilter)

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initRegex = self.inputFilter.getRegex()
        else:
            name = ''
            initRegex = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.regexLineEdit =  QtGui.QLineEdit()
        self.regexLineEdit.setText(initRegex)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('Regex',self.regexLineEdit)


        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            regex = str(self.regexLineEdit.text())
            cRegex = re.compile(regex,re.IGNORECASE)
            returnFilter = DPCLOBFilter(name,self.dpName,self.longDpName,regex)
        except:
            returnFilter = None
        return returnFilter

class DPCLOBFilter(DPFilter):
    def __init__(self,name,dpName,longDpName,regex):
        """Constructor."""

        super(DPCLOBFilter,self).__init__(name)
        self.dpName = dpName
        self.longDpName = longDpName
        self.regex = regex

        self.cRegex = re.compile(self.regex,re.IGNORECASE)



    def getToolTip(self):
        """Gets tooltip."""
        
        return str(self.regex)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches self.regex."""
        
        #get values from xml 
        clobElements = xml.xpath('./'+self.dpName+'/DP_CLOB1')
        #print 'self.value',self.value

        if len(clobElements) > 0:
            #check all values for a match.
            for e in clobElements:
                if self.cRegex.search(e.text):
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False

        else:
            return False
    
    def getRegex(self):
        """Return regex as string."""
        return self.regex

class makeDPVARCHAR2Filter(makeDPCLOBFilter):

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            regex = str(self.regexLineEdit.text())
            cRegex = re.compile(regex,re.IGNORECASE)
            returnFilter = DPVARCHAR2Filter(name,self.dpName,self.longDpName,regex)
        except:
            returnFilter = None
        return returnFilter


class DPVARCHAR2Filter(DPCLOBFilter):
    """VARCHAR2 filter class."""

    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches self.regex."""
        
        #get values from xml 
        varcharElements = xml.xpath('./'+self.dpName+'/DP_QUAL_VALUE1')
        #print 'self.value',self.value

        if len(varcharElements) > 0:
            #check all values for a match.
            for e in varcharElements:
                if self.cRegex.search(e.text):
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False
        else:
            return False
   



class makeDPYNFilter(makeFilter):
    """Makes a Y/N filter."""
    def __init__(self,parent,dpName,longDpName,inputFilter=None):
        """Constructor."""
        super(makeDPYNFilter,self).__init__(parent,inputFilter)
        self.dpName = dpName
        self.longDpName = longDpName

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.valueCombo = QtGui.QComboBox()
        self.valueCombo.addItem('Y','Y')
        self.valueCombo.addItem('N','N')
        self.valueCombo.addItem('Uncertain','Uncertain')

        if initValue == 'Y':
            self.valueCombo.setCurrentIndex(0)
        elif initValue == 'N':
            self.valueCombo.setCurrentIndex(1)
        elif initValue == 'Uncertain':
            self.valueCombo.setCurrentIndex(2)


        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('value',self.valueCombo)


        self.setLayout(flo)



    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            value = str(self.valueCombo.currentText())
            returnFilter = DPYNFilter(name,self.dpName,self.longDpName,value)
        except:
            returnFilter = None
        return returnFilter

class DPYNFilter(DPFilter):
    """DP Y/N filter class."""
    def __init__(self,name,dpName,longDpName,value):
        """Constructor."""

        super(DPYNFilter,self).__init__(name)
        self.dpName = dpName
        self.longDpName = longDpName
        self.value = value

    def getToolTip(self):
        """Returns tooltip."""
        return str(self.value)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches \'value\'"""
        
        #get values from xml 
        qvElements = xml.xpath('./'+self.dpName+'/DP_QUAL_VALUE1')
        #print 'self.value',self.value

        if len(qvElements) > 0:
            #check all values for a match.
            for e in qvElements:
                if e.text == self.value:
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False
        else:
            return False
    
    def getValue(self):
        """Return regex as string."""
        return self.value


class DPVARCHAR2CHOICEFilter(DPYNFilter):
    """DP_VARCHAR2 filter where there are a list of allowed values."""
    def __init__(self,name,dpName,longDpName,value):
        super(DPVARCHAR2CHOICEFilter,self).__init__(name,dpName,longDpName,value)


class makeDPVARCHAR2CHOICEFilter(makeDPVARCHAR2Filter):
    """Makes DPVARCHAR2CHOICE filter."""
    def __init__(self,parent,dpName,longDpName,structureDict,inputFilter=None):
        """Constructor."""
        self.structureDict = structureDict
        super(makeDPVARCHAR2CHOICEFilter,self).__init__(parent,dpName,longDpName,inputFilter)


    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.valueCombo = QtGui.QComboBox()
        self.valueCombo.addItems(\
                list(self.structureDict[self.dpName]['choices']))

        if initValue in list(self.structureDict[self.dpName]['choices']):
            idx = self.valueCombo.findText(initValue,QtCore.Qt.MatchExactly)
            self.valueCombo.setCurrentIndex(idx)
        else:
            self.valueCombo.setCurrentIndex(0)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('value',self.valueCombo)


        self.setLayout(flo)


    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            value = str(self.valueCombo.currentText())
            returnFilter = DPVARCHAR2CHOICEFilter(name,self.dpName,self.longDpName,value)
        except:
            returnFilter = None
        return returnFilter

class DPINTEGERFilter(DPFilter):
    """DPINTEGER filter."""
    def __init__(self,name,dpName,longDpName,integer):
        """Constructor."""

        super(DPINTEGERFilter,self).__init__(name)
        self.dpName = dpName
        self.longDpName = longDpName
        self.integer = integer

    def getToolTip(self):
        """Returns tooltip."""
        
        return str(self.integer)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches self.integer."""
        
        #get values from xml 
        integerElements = xml.xpath('./'+self.dpName+'/DP_INTEGER1')
        #print 'self.value',self.value

        if len(integerElements) > 0:
            #check all values for a match.
            for e in integerElements:
                if str(e.text) == self.integer:
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False

        else:
            return False
    
    def getInteger(self):
        """Return integer as string."""
        return self.integer

class makeDPINTEGERFilter(makeFilter):
    """Makes an integer filter."""
    def __init__(self,parent,dpName,longDpName,inputFilter=None):
        """Constructor."""
        self.dpName = dpName
        self.longDpName = longDpName
        super(makeDPINTEGERFilter,self).__init__(parent,inputFilter)

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initInteger = self.inputFilter.getInteger()
        else:
            name = ''
            initInteger = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.intLineEdit =  QtGui.QLineEdit()
        intV = intBlankValidator()
        self.intLineEdit.setValidator(intV) 
        self.intLineEdit.setText(initInteger)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('Integer',self.intLineEdit)


        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())
            outInteger = str(self.intLineEdit.text())

            returnFilter = DPINTEGERFilter(name,self.dpName,self.longDpName,outInteger)
        except:
            returnFilter = None
        return returnFilter

class makeDPVARCHAR2Filter(makeDPCLOBFilter):
    """Makes DP_VARCHAR2 filter."""

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            regex = str(self.regexLineEdit.text())
            cRegex = re.compile(regex,re.IGNORECASE)
            returnFilter = DPVARCHAR2Filter(name,self.dpName,self.longDpName,regex)
        except:
            returnFilter = None
        return returnFilter


class DPVARCHAR2Filter(DPCLOBFilter):
    """DP_VARCHAR2 filter."""

    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches self.regex."""
        
        #get values from xml 
        varcharElements = xml.xpath('./'+self.dpName+'/DP_QUAL_VALUE1')
        #print 'self.value',self.value

        if len(varcharElements) > 0:
            #check all values for a match.
            for e in varcharElements:
                if self.cRegex.search(e.text):
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False
        else:
            return False
   



class makeDPYNFilter(makeFilter):
    """Makes a Y/N filter."""
    def __init__(self,parent,dpName,longDpName,inputFilter=None):
        """Constructor."""
        super(makeDPYNFilter,self).__init__(parent,inputFilter)
        self.dpName = dpName
        self.longDpName = longDpName

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.valueCombo = QtGui.QComboBox()
        self.valueCombo.addItem('Y','Y')
        self.valueCombo.addItem('N','N')
        self.valueCombo.addItem('Uncertain','Uncertain')

        if initValue == 'Y':
            self.valueCombo.setCurrentIndex(0)
        elif initValue == 'N':
            self.valueCombo.setCurrentIndex(1)
        elif initValue == 'Uncertain':
            self.valueCombo.setCurrentIndex(2)


        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('value',self.valueCombo)


        self.setLayout(flo)



    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            value = str(self.valueCombo.currentText())
            returnFilter = DPYNFilter(name,self.dpName,self.longDpName,value)
        except:
            returnFilter = None
        return returnFilter

class DPYNFilter(DPFilter):
    """DP Y/N filter."""
    def __init__(self,name,dpName,longDpName,value):
        """Constructor."""

        super(DPYNFilter,self).__init__(name)
        self.dpName = dpName
        self.longDpName = longDpName
        self.value = value

    def getToolTip(self):
        """Returns tooltip."""
        return str(self.value)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value that matches \'value\'"""
        
        #get values from xml 
        qvElements = xml.xpath('./'+self.dpName+'/DP_QUAL_VALUE1')
        #print 'self.value',self.value

        if len(qvElements) > 0:
            #check all values for a match.
            for e in qvElements:
                if e.text == self.value:
                    if self.applyQualifierFilters(e.getparent()):
                        return True
            return False
        else:
            return False
    
    def getValue(self):
        """Return regex as string."""
        return self.value


class DPVARCHAR2CHOICEFilter(DPYNFilter):
    """DP_VARCHAR2 filter. Duplicate?"""
    def __init__(self,name,dpName,longDpName,value):
        super(DPVARCHAR2CHOICEFilter,self).__init__(name,dpName,longDpName,value)


class makeDPVARCHAR2CHOICEFilter(makeDPVARCHAR2Filter):
    """Make DP_VARCHAR2choice filter. Duplicate?"""
    def __init__(self,parent,dpName,longDpName,structureDict,inputFilter=None):
        """Constructor."""
        self.structureDict = structureDict
        super(makeDPVARCHAR2CHOICEFilter,self).__init__(parent,dpName,longDpName,inputFilter)


    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.valueCombo = QtGui.QComboBox()
        self.valueCombo.addItems(\
                list(self.structureDict[self.dpName]['choices']))

        if initValue in list(self.structureDict[self.dpName]['choices']):
            idx = self.valueCombo.findText(initValue,QtCore.Qt.MatchExactly)
            self.valueCombo.setCurrentIndex(idx)
        else:
            self.valueCombo.setCurrentIndex(0)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('value',self.valueCombo)


        self.setLayout(flo)


    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            #compile regex to see it's valid.
            value = str(self.valueCombo.currentText())
            returnFilter = DPVARCHAR2CHOICEFilter(name,self.dpName,self.longDpName,value)
        except:
            returnFilter = None
        return returnFilter


class DPFloatFilter(DPFilter):
    """DP_FLOAT filter."""
    def __init__(self,name,dpName,longDpName,dataDict):
        """Constructor."""

        self.dpName = dpName
        self.longDpName = longDpName
        self.dataDict = dataDict
        super(DPFloatFilter,self).__init__(name)

    def getToolTip(self):
        """Returns tooltip."""
        return str(self.dataDict)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value 
        that matches the range given in dataDict."""

        def unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2):
            #Do some initial sanity checking.
            if DP_UNIT_CODE1 != None and DP_UNIT_CODE2 != None:
                if DP_UNIT_CODE1 != DP_UNIT_CODE2:
                    #print 'DP_UNIT_CODE1 '+DP_UNIT_CODE1 +' != DP_UNIT_CODE2 '+DP_UNIT_CODE2 
                    return False
            if self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] != None:
                if self.dataDict['gtUnit'] != self.dataDict['ltUnit']:
                    #print 'self.dataDict[\'gtUnit\'] '+self.dataDict['gtUnit']+\
                    #        ' != self.dataDict[\'ltUnit\'] '+self.dataDict['ltUnit']
                    return False
                
            #Get data unit
            if DP_UNIT_CODE1 != None and DP_UNIT_CODE2 ==None:
                dataUnit = DP_UNIT_CODE1
            elif DP_UNIT_CODE1 == None and DP_UNIT_CODE2 != None:
                dataUnit = DP_UNIT_CODE2
            elif DP_UNIT_CODE1 != None and DP_UNIT_CODE2 != None:
                dataUnit = DP_UNIT_CODE1 # both equal - checked above.
            elif DP_UNIT_CODE1 == None and DP_UNIT_CODE2 == None:
                dataUnit = None
            else:
                print 'unitFilter: something is wrong'
                exit()

            #get dialog unit.

            if self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] == None:
                dialogUnit = self.dataDict['gtUnit']
            elif self.dataDict['gtUnit'] == None and self.dataDict['ltUnit'] != None:
                dialogUnit = self.dataDict['ltUnit']
            elif self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] != None:
                dialogUnit = self.dataDict['ltUnit'] # both equal - checked above.
            elif self.dataDict['gtUnit'] == None and self.dataDict['ltUnit'] == None:
                dialogUnit = None
            else:
                print 'unitFilter: dialog units: something is wrong'
                exit()

            if dataUnit == dialogUnit:
                return True
            else:
                return False


        def gtFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
                DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2):
            """Takes self.dataDict and the xml values and determines whether 
            the values match the filter. """
            
            #print 'gtFilter', DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
            #                        DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2

            
            unitResult = unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2) 
            #print 'unitResult',unitResult
            
            if (self.dataDict['gtCombo'] == '>' or self.dataDict['gtCombo'] \
                    == '>=') and self.dataDict['fuzzy'] == False:
                #If xml data has a lower bound and it's greater than what's in
                #the combo box, return True. Otherwise, return False
                #This works if DP_QUAL_VALUE1 is any one of [None],>,>=
                if DP_QUANT_VALUE1 != None:
                    if self.dataDict['gtCombo'] == '>':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    

class DPFloatFilter(DPFilter):
    """DP_FLOAT filter. Duplicate?"""
    def __init__(self,name,dpName,longDpName,dataDict):
        """Constructor."""

        self.dpName = dpName
        self.longDpName = longDpName
        self.dataDict = dataDict
        super(DPFloatFilter,self).__init__(name)

    def getToolTip(self):
        """Returns tooltip."""
        return str(self.dataDict)


    def applyFilter(self,xml):
        """Checks to see whether any datapoint of type self.dpName has a value 
        that matches the range given in dataDict."""

        def unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2):
            #Do some initial sanity checking.
            if DP_UNIT_CODE1 != None and DP_UNIT_CODE2 != None:
                if DP_UNIT_CODE1 != DP_UNIT_CODE2:
                    #print 'DP_UNIT_CODE1 '+DP_UNIT_CODE1 +' != DP_UNIT_CODE2 '+DP_UNIT_CODE2 
                    return False
            if self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] != None:
                if self.dataDict['gtUnit'] != self.dataDict['ltUnit']:
                    #print 'self.dataDict[\'gtUnit\'] '+self.dataDict['gtUnit']+\
                    #        ' != self.dataDict[\'ltUnit\'] '+self.dataDict['ltUnit']
                    return False
                
            #Get data unit
            if DP_UNIT_CODE1 != None and DP_UNIT_CODE2 ==None:
                dataUnit = DP_UNIT_CODE1
            elif DP_UNIT_CODE1 == None and DP_UNIT_CODE2 != None:
                dataUnit = DP_UNIT_CODE2
            elif DP_UNIT_CODE1 != None and DP_UNIT_CODE2 != None:
                dataUnit = DP_UNIT_CODE1 # both equal - checked above.
            elif DP_UNIT_CODE1 == None and DP_UNIT_CODE2 == None:
                dataUnit = None
            else:
                print 'unitFilter: something is wrong'
                exit()

            #get dialog unit.

            if self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] == None:
                dialogUnit = self.dataDict['gtUnit']
            elif self.dataDict['gtUnit'] == None and self.dataDict['ltUnit'] != None:
                dialogUnit = self.dataDict['ltUnit']
            elif self.dataDict['gtUnit'] != None and self.dataDict['ltUnit'] != None:
                dialogUnit = self.dataDict['ltUnit'] # both equal - checked above.
            elif self.dataDict['gtUnit'] == None and self.dataDict['ltUnit'] == None:
                dialogUnit = None
            else:
                print 'unitFilter: dialog units: something is wrong'
                exit()

            if dataUnit == dialogUnit:
                return True
            else:
                return False


        def gtFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
                DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2):
            """Takes self.dataDict and the xml values and determines whether 
            the values match the filter. """
            
            #print 'gtFilter', DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
            #                        DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2

            
            unitResult = unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2) 
            #print 'unitResult',unitResult
            
            if (self.dataDict['gtCombo'] == '>' or self.dataDict['gtCombo'] \
                    == '>=') and self.dataDict['fuzzy'] == False:
                #If xml data has a lower bound and it's greater than what's in
                #the combo box, return True. Otherwise, return False
                #This works if DP_QUAL_VALUE1 is any one of [None],>,>=
                if DP_QUANT_VALUE1 != None:
                    if self.dataDict['gtCombo'] == '>':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['gtCombo'] == '>=':
                        if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                else:
                    dataResult =  False

            elif (self.dataDict['gtCombo'] == '>' or self.dataDict['gtCombo'] \
                    == '>=') and self.dataDict['fuzzy'] == True:
                #Three cases. Xml data is of type: 
                #1. = (single number)
                #2. < / <=
                #3. > / >=
                #4. Range

                if (DP_QUAL_VALUE1 == None and DP_QUANT_VALUE1 != None \
                        and DP_QUAL_VALUE2 == None and DP_QUANT_VALUE2 == None) or\
                        ((DP_QUAL_VALUE1 == '&gt;' or DP_QUAL_VALUE1 == '&gt;=')\
                        and DP_QUANT_VALUE1 != None and DP_QUAL_VALUE2 == None and \
                        DP_QUANT_VALUE2 == None):
                    #Case 1 or Case 3

                    #Case 1
                    if DP_QUAL_VALUE1 == None:
                        if self.dataDict['gtCombo'] == '>':
                            if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif self.dataDict['gtCombo'] == '>=':
                            if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        else:
                            print 'self.dataDict[\'gtCombo\'] should be > or >=.'
                            exit()
                    else:
                        #Case 3:
                        #Any two ranges where the right bound is infinity
                        #have a region where they overlap.
                        dataResult = True

                elif (DP_QUAL_VALUE1 == None and DP_QUANT_VALUE1 == None and\
                        (DP_QUAL_VALUE2 == '&lt;' or DP_QUAL_VALUE2 == '&lt;=')\
                        and DP_QUANT_VALUE2 != None) or \
                        (DP_QUAL_VALUE1 != None and DP_QUANT_VALUE1 != None\
                        and DP_QUAL_VALUE2 != None and DP_QUANT_VALUE2 != None):
                    #Case 2 or case 4:
                    if self.dataDict['gtCombo'] == '>':
                        if DP_QUANT_VALUE2 > self.dataDict['gtLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['gtCombo'] == '>=':
                        if DP_QUANT_VALUE2 >= self.dataDict['gtLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False

                else:
                    print '???????'
                    exit()
            else:
                print 'dpFloatFilter.gtFilter:???'
                exit()

            return dataResult and unitResult

        def eqFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
                DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2):
            #print 'eqFilter'
            #print 'DP_QUANT_VALUE1',DP_QUANT_VALUE1
            #print 'DP_QUAL_VALUE1',DP_QUAL_VALUE1
            #print 'DP_UNIT_CODE1',DP_UNIT_CODE1
            #print 'DP_QUANT_VALUE2',DP_QUANT_VALUE2
            #print 'DP_QUAL_VALUE2',DP_QUAL_VALUE2
            #print 'DP_UNIT_CODE2',DP_UNIT_CODE2
            #print self.dataDict

            unitResult = unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2) 
            #print 'unitResult',unitResult

            

            if self.dataDict['fuzzy'] == False:
                #print 'fuzzy == False'
                #exact hit.
                if DP_QUANT_VALUE1 == self.dataDict['gtLineEdit'] and \
                        DP_QUAL_VALUE1 == None and DP_QUANT_VALUE2 == None \
                        and DP_QUAL_VALUE2 == None and DP_QUANT_VALUE1 != None:
                    dataResult = True
                else:
                    dataResult = False
            elif self.dataDict['fuzzy'] == True:
                #fuzzy hit.
                #Possibilities in xml:
                #1. single value
                #2. >/>=
                #3. </<=
                #4. Range.

                if DP_QUANT_VALUE1 != None and DP_QUANT_VALUE2 == None\
                        and DP_QUAL_VALUE1 == None and DP_QUAL_VALUE2 == None:
                    #Case 1.

                    if DP_QUANT_VALUE1 == self.dataDict['gtLineEdit']:
                        dataResult = True
                    else:
                        dataResult = False
                elif DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == None:
                    #Case 2.1
                    if DP_QUANT_VALUE1 < self.dataDict['gtLineEdit']:
                        dataResult = True
                    else:
                        dataResult = False
                elif DP_QUAL_VALUE1 == '&gt;=' and DP_QUAL_VALUE2 == None:
                    #Case 2.2
                    if DP_QUANT_VALUE1 <= self.dataDict['gtLineEdit']:
                        dataResult = True
                    else:
                        dataResult = False

                elif DP_QUAL_VALUE1 == None and DP_QUAL_VALUE2 == '&lt;':
                    #Case 3.1
                    if DP_QUANT_VALUE2 > self.dataDict['ltLineEdit']:
                        dataResult = True
                    else:
                        dataResult = False

                elif DP_QUAL_VALUE1 == None and DP_QUAL_VALUE2 == '&lt;=':
                    #Case 3.2
                    if DP_QUANT_VALUE2 >= self.dataDict['ltLineEdit']:
                        dataResult =  True
                    else:
                        dataResult = False

                elif DP_QUAL_VALUE1 != None and DP_QUAL_VALUE2 != None:
                    #Case 4.

                    if DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == '&lt;':
                        if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit'] and\
                                DP_QUANT_VALUE2 > self.DataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult =  False
                    elif DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == '&lt;=':
                        if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit'] and\
                                DP_QUANT_VALUE2 >= self.DataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult =  False
                    elif DP_QUAL_VALUE1 == '&gt;=' and DP_QUAL_VALUE2 == '&lt;':
                        if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit'] and\
                                DP_QUANT_VALUE2 > self.DataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult = False
                    elif DP_QUAL_VALUE1 == '&gt;=' and DP_QUAL_VALUE2 == '&lt;=':
                        if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit'] and\
                                DP_QUANT_VALUE2 >= self.DataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult = False
                    else:
                        print 'Something is terribly wrong.'
                        exit()
                else:
                    print 'Something is wrong.'
                    exit()
            return dataResult and unitResult

        def ltFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
                DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2):
            #print 'ltFilter'
            #print 'DP_QUANT_VALUE1',DP_QUANT_VALUE1
            #print 'DP_QUAL_VALUE1',DP_QUAL_VALUE1
            #print 'DP_UNIT_CODE1',DP_UNIT_CODE1
            #print 'DP_QUANT_VALUE2',DP_QUANT_VALUE2
            #print 'DP_QUAL_VALUE2',DP_QUAL_VALUE2
            #print 'DP_UNIT_CODE2',DP_UNIT_CODE2
            #print self.dataDict
            
            #The following also returns True if both DP_UNIT_CODE1 and 
            #self.dataDict['gtUnit'] are equal to None.
            #if DP_UNIT_CODE1 == self.dataDict['gtUnit'] and \
            #        DP_UNIT_CODE2 == self.dataDict['ltUnit']:
            #    unitResult = True
            #else:
            #    unitResult = False
            #print 'unitResult',unitResult
            unitResult = unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2) 
            #print 'unitResult',unitResult

            if (self.dataDict['ltCombo'] == '<' or self.dataDict['ltCombo'] \
                    == '<=') and self.dataDict['fuzzy'] == False:
                #If xml data has an upper bound and it's less than what's in
                #the combo box, return True. Otherwise, return False
                if DP_QUANT_VALUE2 != None:
                    #There is a right hand bound.
                    if self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE2 < self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE2 <= self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                else:
                    #No right hand bound.
                    if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None:
                        #=
                        #print 'xml equals'
                        if self.dataDict['ltCombo'] == '<':
                            #print DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']
                            if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif self.dataDict['ltCombo'] == '<=':
                            if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        else:
                            print 'self.dataDict[\'ltCombo\'] should be < or <=.'
                            exit()
                    else:
                        dataResult =  False



            elif (self.dataDict['ltCombo'] == '<' or self.dataDict['ltCombo'] \
                    == '<=') and self.dataDict['fuzzy'] == True:

                #Three cases. Xml data is of type: 
                #1. = (single number)
                #2. > / >=
                #3. < / <=
                #4. Range

                if (DP_QUAL_VALUE1 == None and DP_QUANT_VALUE1 != None \
                        and DP_QUAL_VALUE2 == None and DP_QUANT_VALUE2 == None) or\
                        ((DP_QUAL_VALUE2 == '&lt;' or DP_QUAL_VALUE2 == '&lt;=')\
                        and DP_QUANT_VALUE2 != None and DP_QUAL_VALUE1 == None and \
                        DP_QUANT_VALUE1 == None):
                    #Case 1 or Case 3

                    #Case 1
                    if DP_QUAL_VALUE1 == None:
                        if self.dataDict['ltCombo'] == '<':
                            if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif self.dataDict['ltCombo'] == '<=':
                            if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        else:
                            print 'self.dataDict[\'ltCombo\'] should be < or <=.'
                            exit()
                    else:
                        #Case 3:
                        #Any two ranges where the left bound is infinity
                        #have a region where they overlap.
                        dataResult = True

                elif (DP_QUAL_VALUE2 == None and DP_QUANT_VALUE2 == None and\
                        (DP_QUAL_VALUE1 == '&gt;' or DP_QUAL_VALUE1 == '&gt;=')\
                        and DP_QUANT_VALUE1 != None) or \
                        (DP_QUAL_VALUE1 != None and DP_QUANT_VALUE1 != None\
                        and DP_QUAL_VALUE2 != None and DP_QUANT_VALUE2 != None):
                    #Case 2 or case 4:
                    if self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False

                else:
                    print '???????'
                    



            else:
                print 'ltFilter: something went awry.'
                exit()
            
            return dataResult and unitResult
        

        def rangeFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,DP_QUAL_VALUE1,\
                DP_QUANT_VALUE2,DP_UNIT_CODE2,DP_QUAL_VALUE2):
            #print 'rangeFilter'
            #print 'DP_QUANT_VALUE1',DP_QUANT_VALUE1
            #print 'DP_QUAL_VALUE1',DP_QUAL_VALUE1
            #print 'DP_UNIT_CODE1',DP_UNIT_CODE1
            #print 'DP_QUANT_VALUE2',DP_QUANT_VALUE2
            #print 'DP_QUAL_VALUE2',DP_QUAL_VALUE2
            #print 'DP_UNIT_CODE2',DP_UNIT_CODE2
            #print self.dataDict

            unitResult = unitFilter(DP_UNIT_CODE1,DP_UNIT_CODE2) 
            #print 'unitResult',unitResult
            
            if self.dataDict['fuzzy'] == False:
                #Only True if data is a single number in the range
                #Or a range where both the top and bottom numbers are 
                #in the range. False otherwise.

                if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None and \
                        DP_QUANT_VALUE2 == None and DP_QUAL_VALUE1 == None:
                    if self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] and\
                                DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] and\
                                DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] and\
                                DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] and\
                                DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False
                    else:
                        print 'rangeFilter: something is wrong.'
                        exit()
                    

                elif DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 != None and \
                        DP_QUANT_VALUE2 != None and DP_QUAL_VALUE1 != None:
                    if self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE2 > self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE1 < self.dataDict['ltLineEdit'] \
                                and DP_QUANT_VALUE2 < self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult = False

                    elif self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE2 > self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit'] \
                                and DP_QUANT_VALUE2 <= self.dataDict['ltLineEdit']:
                            dataResult = True
                        else:
                            dataResult =  False
                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE2 >= self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE1 < self.dataDict['ltLineEdit'] \
                                and DP_QUANT_VALUE2 < self.dataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult =  False

                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE2 >= self.dataDict['gtLineEdit'] \
                                and DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit'] \
                                and DP_QUANT_VALUE2 <= self.dataDict['ltLineEdit']:
                            dataResult =  True
                        else:
                            dataResult = False

            elif self.dataDict['fuzzy'] == True:
                    if self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None\
                                and DP_QUANT_VALUE2 == None and \
                                DP_QUAL_VALUE2 == None:
                            # =
                            if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] \
                                    and DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                 or DP_QUAL_VALUE1 == '&gt;=') and \
                                 DP_QUANT_VALUE2 == None and \
                                 DP_QUAL_VALUE2 == None:
                            # > / >=
                            if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False

                        elif DP_QUANT_VALUE1 == None and DP_QUAL_VALUE1 == None and\
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            # < / <=
                            if DP_QUANT_VALUE2 > self.dataDict['gtLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                or DP_QUAL_VALUE1 == '&gt;=') and \
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            #Range
                            if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                    DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                dataResult = False
                            else:
                                dataResult = True
                            
                    elif self.dataDict['gtCombo'] == '>' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None\
                                and DP_QUANT_VALUE2 == None and \
                                DP_QUAL_VALUE2 == None:
                            # =
                            if DP_QUANT_VALUE1 > self.dataDict['gtLineEdit'] \
                                    and DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                 or DP_QUAL_VALUE1 == '&gt;=') and \
                                 DP_QUANT_VALUE2 == None and \
                                 DP_QUAL_VALUE2 == None:
                            if DP_QUAL_VALUE1 == '&gt;':
                                # > 
                                if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False
                            elif DP_QUAL_VALUE1 == '&gt;=':
                                # >=
                                if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False

                        elif DP_QUANT_VALUE1 == None and DP_QUAL_VALUE1 == None and\
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            # < / <=
                            if DP_QUANT_VALUE2 > self.dataDict['gtLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                or DP_QUAL_VALUE1 == '&gt;=') and \
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            #Range
                            if DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True
                            elif DP_QUAL_VALUE2 == '&lt;=':
                                if DP_QUANT_VALUE2 < self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True

                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<':
                        if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None\
                                and DP_QUANT_VALUE2 == None and \
                                DP_QUAL_VALUE2 == None:
                            # =
                            if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] \
                                    and DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                 or DP_QUAL_VALUE1 == '&gt;=') and \
                                 DP_QUANT_VALUE2 == None and \
                                 DP_QUAL_VALUE2 == None:
                            # > / >=
                            if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False

                        elif DP_QUANT_VALUE1 == None and DP_QUAL_VALUE1 == None and\
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            # < / <=
                            if DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 > self.dataDict['gtLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False
                            elif DP_QUAL_VALUE2 == '&lt;=':
                                if DP_QUANT_VALUE2 >= self.dataDict['gtLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False


                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                or DP_QUAL_VALUE1 == '&gt;=') and \
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            #Range
                            if DP_QUAL_VALUE1 == '&gt;':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True
                            elif DP_QUAL_VALUE1 == '&gt;=':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 > self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True



                    elif self.dataDict['gtCombo'] == '>=' \
                            and self.dataDict['ltCombo'] == '<=':
                        if DP_QUANT_VALUE1 != None and DP_QUAL_VALUE1 == None\
                                and DP_QUANT_VALUE2 == None and \
                                DP_QUAL_VALUE2 == None:
                            # =
                            if DP_QUANT_VALUE1 >= self.dataDict['gtLineEdit'] \
                                    and DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                dataResult = True
                            else:
                                dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                 or DP_QUAL_VALUE1 == '&gt;=') and \
                                 DP_QUANT_VALUE2 == None and \
                                 DP_QUAL_VALUE2 == None:
                            if DP_QUAL_VALUE1 == '&gt;':
                                # > 
                                if DP_QUANT_VALUE1 < self.dataDict['ltLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False
                            elif DP_QUAL_VALUE1 == '&gt;=':
                                # >=
                                if DP_QUANT_VALUE1 <= self.dataDict['ltLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False




                        elif DP_QUANT_VALUE1 == None and DP_QUAL_VALUE1 == None and\
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            # < / <=
                            if DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 > self.dataDict['gtLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False
                            elif DP_QUAL_VALUE2 == '&lt;=':
                                if DP_QUANT_VALUE2 >= self.dataDict['gtLineEdit']:
                                    dataResult = True
                                else:
                                    dataResult = False
                        elif DP_QUANT_VALUE1 != None and (DP_QUAL_VALUE1 == '&gt;'\
                                or DP_QUAL_VALUE1 == '&gt;=') and \
                                DP_QUANT_VALUE2 != None and (DP_QUAL_VALUE2 =='&lt;'\
                                or DP_QUAL_VALUE2 == '&lt;='):
                            #Range
                            if DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True
                            elif DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == '&lt;=':
                                if DP_QUANT_VALUE2 < self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True

                            elif DP_QUAL_VALUE1 == '&gt;' and DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 >= self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True
                            elif DP_QUAL_VALUE1 == '&gt;=' and DP_QUAL_VALUE2 == '&lt;':
                                if DP_QUANT_VALUE2 <= self.dataDict['gtLineEdit'] or\
                                        DP_QUANT_VALUE1 > self.dataDict['ltLineEdit']:
                                    dataResult = False
                                else:
                                    dataResult = True


            return dataResult and unitResult
            

        
        #print 'DPFloatFilter'
        #get applicable xml elements
        dpElements = xml.xpath('./'+self.dpName)

        if len(dpElements) > 0:
            for e in dpElements:
                #check datapoint for a match.
                DP_QUANT_VALUE1 = e.xpath('./DP_QUANT_VALUE1')
                DP_UNIT_CODE1 = e.xpath('./DP_UNIT_CODE1')
                DP_QUAL_VALUE1 = e.xpath('./DP_QUAL_VALUE1')
                DP_QUANT_VALUE2 = e.xpath('./DP_QUANT_VALUE2')
                DP_UNIT_CODE2 = e.xpath('./DP_UNIT_CODE2')
                DP_QUAL_VALUE2 = e.xpath('./DP_QUAL_VALUE2')

                if len(DP_QUANT_VALUE1) > 1 or len(DP_UNIT_CODE1) > 1 or len(DP_QUANT_VALUE2) > 1 or len(DP_UNIT_CODE2) > 1 or len(DP_QUAL_VALUE1) > 1 or len(DP_QUAL_VALUE2) > 1:
                    print 'DPFloatFilter: ???'
                    exit()

                if len(DP_QUANT_VALUE1) == 1:
                    DP_QUANT_VALUE1 = float(DP_QUANT_VALUE1[0].text)
                else:
                    #len(DP_QUANT_VALUE1) == 0
                    DP_QUANT_VALUE1 = None

                if len(DP_UNIT_CODE1) == 1:
                    DP_UNIT_CODE1 = DP_UNIT_CODE1[0].text
                else:
                    DP_UNIT_CODE1 = None
                
                if len(DP_QUAL_VALUE1) == 1:
                    DP_QUAL_VALUE1 = DP_QUAL_VALUE1[0].text
                else:
                    DP_QUAL_VALUE1 = None

                if len(DP_QUANT_VALUE2) == 1:
                    DP_QUANT_VALUE2 = float(DP_QUANT_VALUE2[0].text)
                else:
                    DP_QUANT_VALUE2 = None

                if len(DP_UNIT_CODE2) == 1:
                    DP_UNIT_CODE2 = DP_UNIT_CODE2[0].text
                else:
                    DP_UNIT_CODE2 = None
                
                if len(DP_QUAL_VALUE2) == 1:
                    DP_QUAL_VALUE2 = DP_QUAL_VALUE2[0].text
                else:
                    DP_QUAL_VALUE2 = None


                hit = None
                if (self.dataDict['gtCombo'] == '>' or \
                        self.dataDict['gtCombo'] == '>=') and \
                        self.dataDict['ltCombo'] == None:
                    #print 'running gtFilter'
                    hit = gtFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,\
                        DP_QUAL_VALUE1,DP_QUANT_VALUE2,DP_UNIT_CODE2,\
                        DP_QUAL_VALUE2)
                elif self.dataDict['gtCombo'] == '=' and \
                        self.dataDict['ltCombo']== None:
                    #print 'running eqFilter'
                    hit = eqFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,\
                        DP_QUAL_VALUE1,DP_QUANT_VALUE2,DP_UNIT_CODE2,\
                        DP_QUAL_VALUE2)

                elif (self.dataDict['ltCombo'] == '<' or \
                        self.dataDict['ltCombo'] == '<=') and \
                        self.dataDict['gtCombo'] == None:
                    #print 'running ltFilter'
                    hit = ltFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,\
                        DP_QUAL_VALUE1,DP_QUANT_VALUE2,DP_UNIT_CODE2,\
                        DP_QUAL_VALUE2)
                elif (self.dataDict['ltCombo'] == '<' or \
                        self.dataDict['ltCombo'] == '<=') and \
                        (self.dataDict['gtCombo'] == '>' or \
                        self.dataDict['gtCombo'] == '>='):
                    hit = rangeFilter(DP_QUANT_VALUE1,DP_UNIT_CODE1,\
                        DP_QUAL_VALUE1,DP_QUANT_VALUE2,DP_UNIT_CODE2,\
                        DP_QUAL_VALUE2)


                else:
                    print 'Not implemented yet!'
                    print 'self.dataDict: ',self.dataDict
                    exit()
                    hit = True


                if hit == True:
                    if self.applyQualifierFilters(e):
                        return True
            return False
        else:
            return False
    
    def getValue(self):
        """Return self.dataDict."""
        return self.dataDict





class makeDPFloatFilter(makeFilter):
    """Makes a DPFloat filter."""
    def __init__(self,parent,dpName,longDpName,inputFilter=None):
        """Constructor."""
        self.dpName = dpName
        self.longDpName = longDpName
        super(makeDPFloatFilter,self).__init__(parent,inputFilter)

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = {'ltCombo': None, 'ltLineEdit': None, 'ltUnit': None, 'gtCombo': None, 'gtUnit': None, 'fuzzy': False, 'gtLineEdit': None}

        dv = QtGui.QDoubleValidator()
        
        self.nameLineEdit = QtGui.QLineEdit(self)
        self.nameLineEdit.setText(name)

        self.gtCombo = QtGui.QComboBox()
        self.gtCombo.addItem('>','>')
        self.gtCombo.addItem('>=','>=')
        self.gtCombo.addItem('=','=')
        self.gtLineEdit = QtGui.QLineEdit()
        self.gtLineEdit.setValidator(dv)

        self.ltCombo = QtGui.QComboBox()
        self.ltCombo.addItem('<','<')
        self.ltCombo.addItem('<=','<=')
        self.ltLineEdit = QtGui.QLineEdit(self)
        self.ltLineEdit.setValidator(dv)

        if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
            self.ltUnitCombo = QtGui.QComboBox()
            self.gtUnitCombo = QtGui.QComboBox()
            self.ltUnitCombo.addItems(list(self.parentWidget().\
                    structureDict[self.dpName]['dpUnits']))
            self.gtUnitCombo.addItems(list(self.parentWidget().\
                    structureDict[self.dpName]['dpUnits']))


        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('>/>=/=',self.gtCombo)
        flo.addRow('value',self.gtLineEdit)
        if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
            flo.addRow('Unit',self.gtUnitCombo)
            #make sure we're always looking for the same unit.
            self.ltUnitCombo.currentIndexChanged.connect(self.onDPUnitCode1Changed)
            self.gtUnitCombo.currentIndexChanged.connect(self.onDPUnitCode2Changed)
        self.gtCombo.currentIndexChanged.connect(self.onGtComboChanged)


        flo.addRow('</<=',self.ltCombo)
        flo.addRow('value',self.ltLineEdit)
        if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
            flo.addRow('Unit',self.ltUnitCombo)

        self.fuzzy = QtGui.QCheckBox()
        flo.addRow('Include fuzzy hits',self.fuzzy)

        #populate data from self.inputFilter
        if self.inputFilter != None:
            print initValue
            if initValue['ltLineEdit'] != None:
                self.ltLineEdit.setText(str(initValue['ltLineEdit']))
            if initValue['gtLineEdit'] != None:
                self.gtLineEdit.setText(str(initValue['gtLineEdit']))
            if initValue['ltCombo'] != None:
                ltComboIdx = self.ltCombo.findData(initValue['ltCombo'])
                self.ltCombo.setCurrentIndex(ltComboIdx)
            if initValue['gtCombo'] != None:
                gtComboIdx = self.gtCombo.findData(initValue['gtCombo'])
                self.gtCombo.setCurrentIndex(gtComboIdx)
            if hasattr(self,'ltUnitCombo'):
                if initValue['ltUnit'] != None:
                    ltUnitComboIdx = self.ltUnitCombo.findData(initValue['ltUnit'])
                    self.ltUnitCombo.setCurrentIndex(ltUnitComboIdx)
            if hasattr(self,'gtUnitCombo'):
                if initValue['gtUnit'] != None:
                    gtUnitComboIdx = self.gtUnitCombo.findData(initValue['gtUnit'])
                    self.gtUnitCombo.setCurrentIndex(gtUnitComboIdx)


        self.onGtComboChanged()

        self.setLayout(flo)

    def onDPUnitCode1Changed(self):
        """Sync UnitCode1 and UnitCode2"""
        currentText = self.ltUnitCombo.currentText()
        idx = self.gtUnitCombo.findText(currentText,QtCore.Qt.MatchExactly)
        self.gtUnitCombo.setCurrentIndex(idx)
        

    
    def onDPUnitCode2Changed(self):
        """Sync UnitCode1 and UnitCode2"""
        currentText = self.gtUnitCombo.currentText()
        idx = self.ltUnitCombo.findText(currentText,QtCore.Qt.MatchExactly)
        self.ltUnitCombo.setCurrentIndex(idx)

    def onGtComboChanged(self):
        """Grey out gtUnitCombo, gtLineEdit and gtUnitCombo if ltUnitCombo
        is set to "=". 
        
        Activate them when ltUnitCombo is changed from "=" to
        something else."""

        #print 'onLtComboChanged'
        currentText = str(self.gtCombo.currentText())
        if currentText == '=':
            self.ltCombo.setDisabled(True)
            self.ltLineEdit.setDisabled(True)
            self.ltLineEdit.setText('')
            if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
                self.ltUnitCombo.setDisabled(True)
        else:
            self.ltCombo.setDisabled(False)
            self.ltLineEdit.setDisabled(False)
            if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
                self.ltUnitCombo.setDisabled(False)

                
            

    def makeFilter(self):
        """Makes the filter."""
        name = str(self.nameLineEdit.text())

        #get information from combo boxes etc and put into a dict
        #so that it can be passed to DPFloatFilter
        dataDict = {'ltCombo' : None, 'ltLineEdit' : None, 'ltUnit' : None,\
                'gtCombo' : None, 'gtLineEdit' : None, 'gtUnit' : None, \
                'fuzzy' : None}

        if str(self.ltLineEdit.text()) != '':
            dataDict['ltCombo'] = str(self.ltCombo.currentText())
            dataDict['ltLineEdit'] = float(str(self.ltLineEdit.text()))
            if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
                dataDict['ltUnit'] = str(self.ltUnitCombo.currentText())

        if str(self.gtLineEdit.text()) != '':
            dataDict['gtCombo'] = str(self.gtCombo.currentText())
            dataDict['gtLineEdit'] = float(str(self.gtLineEdit.text()))
            if len(self.parentWidget().structureDict[self.dpName]['dpUnits']) > 0:
                dataDict['gtUnit'] = str(self.gtUnitCombo.currentText())

        
        dataDict['fuzzy'] = self.fuzzy.isChecked()
        
        returnFilter = DPFloatFilter(name,self.dpName,self.longDpName,dataDict)

        return returnFilter



class makeQYNFilter(makeFilter):
    """Makes a Y/N filter."""
    def __init__(self,parent,qName,longQName,inputFilter=None):
        """Constructor."""
        self.qName = qName
        self.longQName = longQName
        super(makeQYNFilter,self).__init__(parent,inputFilter)
        #print 'makeQYNFilter.qName: ',self.qName
        #print 'makeQYNFilter.longQName: ',self.longQName

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initValue = self.inputFilter.getValue()
        else:
            name = ''
            initValue = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.valueCombo = QtGui.QComboBox()
        self.valueCombo.addItem('Y','Y')
        self.valueCombo.addItem('N','N')
        self.valueCombo.addItem('Uncertain','Uncertain')

        if initValue == 'Y':
            self.valueCombo.setCurrentIndex(0)
        elif initValue == 'N':
            self.valueCombo.setCurrentIndex(1)
        elif initValue == 'Uncertain':
            self.valueCombo.setCurrentIndex(2)

        self.mostRecentOnly = QtGui.QCheckBox('Most recent only')

        #set self.mostRecentOnly checkbox to a sensible
        #default. It should be checked for ACTIVE and PIVOTAL
        #endpoints, and unchecked otherwise.

        #The user shouldn't have to change this
        #except in very special circumstances.

        if self.qName == 'ACTIVE' or self.qName == 'PIVOTAL':
            self.mostRecentOnly.setChecked(True)
        else:
            self.mostRecentOnly.setChecked(False)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('value',self.valueCombo)
        flo.addRow('',self.mostRecentOnly)

        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        try:
            name = str(self.nameLineEdit.text())

            value = str(self.valueCombo.currentText())

            mostRecent = self.mostRecentOnly.isChecked()

            returnFilter = QYNFilter(name,self.qName,mostRecent,value)
        except:
            returnFilter = None
        return returnFilter



class QYNFilter(filter):
    """Qualifier Y/N filter."""
    def __init__(self,name,qualifierName,mostRecent,inputValue=None):
        """Constructor.
        
        Note, mostRecent is used to restrict the search to the most recent
        qualifier of the given type.
        
        Note: this will also work for YNPStar-type qualifiers."""

        super(QYNFilter,self).__init__(name)
        self.qualifierName = qualifierName
        #mostRecent must be either True or False
        self.mostRecent = mostRecent
        self.value = inputValue

    def getToolTip(self):
        """Returns tooltip."""
        
        return str(self.value)


    def applyFilter(self,xml):
        """Checks to see whether any qualifier of type self.qualifierName has a value equal to self.value."""
        
        #get values from xml 
        values = [x.text for x in xml.xpath('./'+self.qualifierName+'/Q_QUAL_VALUE1')]
        #print 'values',values
        #print 'self.value',self.value
        #print etree.tostring(xml,pretty_print=True)

        if len(values) > 0:
            if self.mostRecent == True:
                #check most recent value only.
                if values[-1] == self.value:
                    return True
                else:
                    return False
            else:
                #check all values for a match.
                for v in values:
                    if v == self.value:
                        return True
                return False

        return False




class QCLOBFilter(filter):
    """Qualifier CLOB filter."""
    def __init__(self,name,qualifierName,mostRecent,regexp=None,):
        """Constructor.
        
        Note, mostRecent is used to restrict the search to the most recent
        qualifier of the given type."""
        

        super(QCLOBFilter,self).__init__(name)
        self.qualifierName = qualifierName
        #mostRecent must be either True or False
        self.mostRecent = mostRecent

        cs = QtCore.Qt.CaseInsensitive
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        self.regexp=QtCore.QRegExp(regexp,cs,syntax)



    def getToolTip(self):
        """Returns tooltip."""
        
        return str(self.value)


    def applyFilter(self,xml):
        """Checks to see whether any qualifier of type self.qualifierName has a value that matches self.value."""
        
        #get values from xml 
        clobElements = xml.xpath('./'+self.qualifierName+'/Q_CLOB1')

        if len(clobElements) > 0:
            if mostRecent == True:
                #check most recent value only.
                if self.regexp.indexIn(clobElements[-1].text) >= 0:
                    return True
                else:
                    return False
            else:
                #check all values for a match.
                for c in clobElements:
                    if self.regexp.indexIn(c.text) >= 0:
                        return True
                return False

        else:
            return False


class makeQCLOBFilter(makeFilter):
    """Makes a filter to search QCLOB filters ."""
    def __init__(self,parent,qName,longQName,inputFilter=None):
        """Constructor."""
        self.qName = qName
        self.longQName = longQName
        super(makeQCLOBFilter,self).__init__(parent,inputFilter)

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initRegex = self.inputFilter.getRegex()
        else:
            name = ''
            initRegex = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.regexLineEdit =  QtGui.QLineEdit()
        self.regexLineEdit.setText(initRegex)

        self.mostRecentOnly = QtGui.QCheckBox('Most recent only')
        self.mostRecentOnly.setChecked(False)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('Regex',self.regexLineEdit)



        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        #try:
        name = str(self.nameLineEdit.text())

        #compile regex to see it's valid.
        regex = str(self.regexLineEdit.text())
        cRegex = re.compile(regex,re.IGNORECASE)
        mostRecent = self.mostRecentOnly.isChecked()
        returnFilter = QCLOBFilter(name,self.qName,mostRecent,regex)
        #except:
        #    returnFilter = None
        #return returnFilter

        return returnFilter


class QVARCHAR2Filter(filter):
    """Qualifier VARCHAR2 filter."""
    def __init__(self,name,qualifierName,mostRecent,regexp=None,):
        """Constructor.
        
        Note, mostRecent is used to restrict the search to the most recent
        qualifier of the given type."""
        

        super(QVARCHAR2Filter,self).__init__(name)
        self.qualifierName = qualifierName
        #mostRecent must be either True or False
        self.mostRecent = mostRecent

        cs = QtCore.Qt.CaseInsensitive
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        self.regexp=QtCore.QRegExp(regexp,cs,syntax)



    def getToolTip(self):
        """Returns tooltip."""
        
        return str(self.value)


    def applyFilter(self,xml):
        """Checks to see whether any qualifier of type self.qualifierName has a value that matches self.value."""
        
        #get values from xml 
        clobElements = xml.xpath('./'+self.qualifierName+'/Q_QUAL_VALUE1')

        if len(clobElements) > 0:
            if mostRecent == True:
                #check most recent value only.
                if self.regexp.indexIn(clobElements[-1].text) >= 0:
                    return True
                else:
                    return False
            else:
                #check all values for a match.
                for c in clobElements:
                    if self.regexp.indexIn(c.text) >= 0:
                        return True
                return False

        else:
            return False



class makeQVARCHAR2Filter(makeFilter):
    """Makes a filter to search QVARCHAR2 filters ."""
    def __init__(self,parent,qName,longQName,inputFilter=None):
        """Constructor."""
        self.qName = qName
        self.longQName = longQName
        super(makeQVARCHAR2Filter,self).__init__(parent,inputFilter)

    def initUI(self):
        """Initialise user interface."""

        if self.inputFilter != None:
            name = self.inputFilter.getName()
            initRegex = self.inputFilter.getRegex()
        else:
            name = ''
            initRegex = ''

        
        self.nameLineEdit = QtGui.QLineEdit()
        self.nameLineEdit.setText(name)
        self.regexLineEdit =  QtGui.QLineEdit()
        self.regexLineEdit.setText(initRegex)

        self.mostRecentOnly = QtGui.QCheckBox('Most recent only')
        self.mostRecentOnly.setChecked(False)

        flo = QtGui.QFormLayout()
        flo.addRow('Filter name',self.nameLineEdit)
        flo.addRow('Regex',self.regexLineEdit)



        self.setLayout(flo)

    def makeFilter(self):
        """Makes the filter."""
        #try:
        name = str(self.nameLineEdit.text())

        #compile regex to see it's valid.
        regex = str(self.regexLineEdit.text())
        cRegex = re.compile(regex,re.IGNORECASE)
        mostRecent = self.mostRecentOnly.isChecked()
        returnFilter = QVARCHAR2Filter(name,self.qName,mostRecent,regex)
        #except:
        #    returnFilter = None
        #return returnFilter

        return returnFilter

