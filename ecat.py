#!/usr/bin/python
##
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
##
## This file incorporates work covered by the following copyright and permission 
## notice
##
#############################################################################
##
## Copyright (C) 2010 Riverbank Computing Limited.
## Copyright (C) 2010 Nokia Corporation and/or its subsidiary(-ies).
## All rights reserved.
##
## This file is part of the examples of PyQt.
##
## $QT_BEGIN_LICENSE:BSD$
## You may use this file under the terms of the BSD license as follows:
##
## "Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
##     the names of its contributors may be used to endorse or promote
##     products derived from this software without specific prior written
##     permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
## $QT_END_LICENSE$
##
#############################################################################

"""\mainpage ECAT documentation

\section Introduction
The Electronic Chemical Assessment System is designed to allow the user to view 
and add to an xml version of a chemical database. The structure of the xml database
is defined in a relaxNG schema, and when new data are added, the viewer 
queries this schema to determine the allowed data types.

The PyQt model/view framework is used throughout. A PyQt model is generated,
which acts as a wrapper for the underlying xml document. This xml document is a
tree with three levels: chemical elements, datapoint elements and qualifier
elements.  The chemical elements are essentially containers for datapoint
elements. Datapoint elements contain all of the data for a chemical, including
CAS-RN's, physico-chemical and toxicological data, and chunks of text
corresponding to sections of a report. Qualifiers contain metadata for each
datapoint, such as tags for searching, citations and comments.

Three PyQt QListViews display the chemicals, datapoints and qualifiers. Between
each view and the model is a QSortFilterProxyModel, allowing a particular
chemical, datapoint or qualifier or groups thereof to be selected with a few
keystrokes.

The structures of each chemical are displayed in the chemical view where
available. If a chemical has a datapoint of type "SMILES", openbabel is used to
convert the SMILES string to an svg graphic, and this is displayed in the
chemical view.
"""

from PyQt4 import QtCore, QtGui, QtSvg
from lxml import etree
from lxml import objectify
from lxml.html.clean import clean_html
import sys
import os
from pprint import pprint 
import copy
import openbabel
import re
import organiseSVG
import pybel
from search import *
import bz2
from shutil import copyfileobj
import shutil
import gc

import autosave
import time
#import objgraph

#required for tree model to avoid cyclic references and a serious memory leak.
import weakref
import cPickle

try:
    import dbXml
except:
    print 'Couldn\'t import dbXml.'

#This imports the plugins and registers them.
from plugins import *


class intBlankValidator(QtGui.QIntValidator):
    """Integer validator that also allows blanks"""
    def __init__(self):
        super(intBlankValidator,self).__init__()
    
    def validate(self,inString,pos):
        """Returns True if inString matches ^[0-9]+$ or is a blank string, False otherwise."""
        if str(inString) == '':
            return (QtGui.QValidator.Acceptable,pos)
        return super(intBlankValidator,self).validate(inString,pos)


class MyDataPointDelegate(QtGui.QStyledItemDelegate):
    """Delegate to display datapoints."""
    def __init__(self):

        super(MyDataPointDelegate,self).__init__()
    
    def paint(self,painter,option,index):
        """Paints delegate in view."""
        #options = QtGui.QStyleOptionViewItemV4()
    #self.initStyleOption(options,index)
        #print 'paint'

        #if QtCore.QString(u'DP_L1_KINGDOM_CODE') in \
        #        index.data(QtCore.Qt.UserRole).toMap():

        painter.save()

        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))

        c1 = None
        c2 = None
        dd = index.data(QtCore.Qt.UserRole).toPyObject()
        if QtCore.QString(u'DP_CONFIDENTIAL1_YN') in dd.keys():
            c1 = index.data(QtCore.Qt.UserRole).toPyObject()[QtCore.QString(u'DP_CONFIDENTIAL1_YN')]
        if QtCore.QString(u'DP_CONFIDENTIAL2_YN') in dd.keys():
            c2 = index.data(QtCore.Qt.UserRole).toPyObject()[QtCore.QString(u'DP_CONFIDENTIAL2_YN')]

        if c1 == QtCore.QString(u'Y') or c2 == QtCore.QString(u'Y'):
            isConfidential = True
        else:
            isConfidential = False
        

        if option.state & QtGui.QStyle.State_Selected:
            if isConfidential == True:
                #dark red
                brush = QtGui.QBrush(QtGui.QColor('#ff0000'))
                painter.setBrush(brush)
            else:
                #green
                brush = QtGui.QBrush(QtGui.QColor('#66ff71'))
                painter.setBrush(brush)

        else:
            if isConfidential == True:
                #light red
                brush = QtGui.QBrush(QtGui.QColor('#ffb3b3'))
                painter.setBrush(brush)
            else:
                pass
            #print index.data(QtCore.Qt.UserRole).toPyObject()


        painter.drawRect(option.rect)
        painter.setPen(QtGui.QPen(QtCore.Qt.blue))
        dispValue = index.data(QtCore.Qt.DisplayRole)

        if dispValue.isValid():
            align = QtCore.Qt.AlignLeft
            
        painter.translate(option.rect.left(),option.rect.top())
        valueDoc = QtGui.QTextDocument()
        valueDoc.setTextWidth(option.rect.width())
        #valueDoc.setHtml(dispValue+' ('+str(valueDoc.pageSize().width())+','+str(valueDoc.pageSize().height())+','+str(valueDoc.size().width())+','+str(valueDoc.size().height())+')')
        valueDoc.setHtml(dispValue.toPyObject())
        #valueDoc.setHtml(dispValue+' ('+str(valueDoc.size().width())+')')
        valueDoc.drawContents(painter)
            #index.model().emit(QtCore.SIGNAL('layoutChanged'))
            #index.model().emit(QtCore.SIGNAL('dataChanged'))


        painter.restore()

    def sizeHint(self,option,index):
        """Generates sizehint for the view."""
        #print 'sizehint'

        dispValue = index.data(QtCore.Qt.DisplayRole)
    
        if dispValue.isValid():
            align = QtCore.Qt.AlignLeft
            valueDoc = QtGui.QTextDocument()
            valueDoc.setTextWidth(option.rect.width())
            valueDoc.setHtml(dispValue.toPyObject())

            return QtCore.QSize(option.rect.width(),valueDoc.size().height())
        else:
            return  QtCore.QSize(option.rect.width(),200)

class MyQualifierDelegate(QtGui.QStyledItemDelegate):
    """Delegate to display qualifiers."""
    def __init__(self):

        super(MyQualifierDelegate,self).__init__()
    
    def paint(self,painter,option,index):
        """Paints delegate in view."""
        #options = QtGui.QStyleOptionViewItemV4()
    #self.initStyleOption(options,index)
        #print 'paint'
        #if QtCore.QString(u'Q_L1_KINGDOM_CODE') in \
        #        index.data(QtCore.Qt.UserRole).toMap():

        painter.save()

        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))



        c1 = None
        c2 = None
        dd = index.data(QtCore.Qt.UserRole).toPyObject()
        if QtCore.QString(u'DP_CONFIDENTIAL1_YN') in dd.keys():
            c1 = index.data(QtCore.Qt.UserRole).toPyObject()[QtCore.QString(u'DP_CONFIDENTIAL1_YN')]
        if QtCore.QString(u'DP_CONFIDENTIAL2_YN') in dd.keys():
            c2 = index.data(QtCore.Qt.UserRole).toPyObject()[QtCore.QString(u'DP_CONFIDENTIAL2_YN')]

        if c1 == QtCore.QString(u'Y') or c2 == QtCore.QString(u'Y'):
            isConfidential = True
        else:
            isConfidential = False



        if option.state & QtGui.QStyle.State_Selected:
            if isConfidential == True:
                #dark red
                brush = QtGui.QBrush(QtGui.QColor('#ff0000'))
                painter.setBrush(brush)
            else:
                #green
                brush = QtGui.QBrush(QtGui.QColor('#66ff71'))
                painter.setBrush(brush)

        else:
            if isConfidential == True:
                #light red
                brush = QtGui.QBrush(QtGui.QColor('#ffb3b3'))
                painter.setBrush(brush)
            else:
                pass
            

        painter.drawRect(option.rect)
        painter.setPen(QtGui.QPen(QtCore.Qt.blue))
        dispValue = index.data(QtCore.Qt.DisplayRole)

        if dispValue.isValid():
            align = QtCore.Qt.AlignLeft
            
        painter.translate(option.rect.left(),option.rect.top())
        valueDoc = QtGui.QTextDocument()
        valueDoc.setTextWidth(option.rect.width())
        #valueDoc.setHtml(dispValue+' ('+str(valueDoc.pageSize().width())+','+str(valueDoc.pageSize().height())+','+str(valueDoc.size().width())+','+str(valueDoc.size().height())+')')
        valueDoc.setHtml(dispValue.toPyObject())
        #valueDoc.setHtml(dispValue+' ('+str(valueDoc.size().width())+')')
        valueDoc.drawContents(painter)
            #index.model().emit(QtCore.SIGNAL('layoutChanged'))
            #index.model().emit(QtCore.SIGNAL('dataChanged'))


        painter.restore()

    def sizeHint(self,option,index):
        """Generates sizehint for the view."""
        #print 'sizehint'

        dispValue = index.data(QtCore.Qt.DisplayRole)
    
        if dispValue.isValid():
            align = QtCore.Qt.AlignLeft
            valueDoc = QtGui.QTextDocument()
            valueDoc.setTextWidth(option.rect.width())
            valueDoc.setHtml(dispValue.toPyObject())

            return QtCore.QSize(option.rect.width(),valueDoc.size().height())
        else:
            return  QtCore.QSize(option.rect.width(),200)

class MyChemicalDelegate(QtGui.QStyledItemDelegate):
    """Delegate to display chemicals, with svg figure and CAS-RN."""
    def __init__(self):

        self.casrnXPath = etree.XPath('./CASNUMBER[ACTIVE[last()]/Q_QUAL_VALUE1=\'Y\']/DP_INTEGER1')
        #self.smilesXPath = etree.XPath('./SMILES/DP_CLOB1')
        super(MyChemicalDelegate,self).__init__()
    
    def paint(self,painter,option,index):
        """Paints delegate in view."""
        #options = QtGui.QStyleOptionViewItemV4()
    #self.initStyleOption(options,index)

        painter.save()

        painter.setPen(QtGui.QPen(QtCore.Qt.NoPen))
    
        if option.state & QtGui.QStyle.State_Selected:
            brush = QtGui.QBrush(QtGui.QColor('#66ff71'))
            painter.setBrush(brush)

        #get casrns.
        chemical =  index.model().mapToSource(index).internalPointer().element
        casrns = [x.text for x in self.casrnXPath(chemical) ]

        
        for idx in range(len(casrns)):
            #print 'casrn',type(casrns[idx])
            if self.cascheck(casrns[idx]) == False:
                casrns[idx] = '<span style="color:red">'+casrns[idx]+"</span>"
    
        casrnstring = '<b>' + ', '.join(casrns) + '</b>'

        painter.drawRect(option.rect)
        painter.setPen(QtGui.QPen(QtCore.Qt.blue))
        dispValue = casrnstring
    
        if index.isValid():
            align = QtCore.Qt.AlignLeft
            
        painter.translate(option.rect.left(),option.rect.top())
        valueDoc = QtGui.QTextDocument()
        valueDoc.setTextWidth(option.rect.width())
        valueDoc.setHtml(dispValue)

        if index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('svgString')].toPyObject() != '':
            valueDoc = self.addSvg(index,valueDoc)
        valueDoc.drawContents(painter)

        painter.restore()

    def cascheck(self,casrn):
        casrn = casrn.strip()
        # check that basic format is valid
        if not re.search('^[0-9]{5}[0-9]*$',casrn) :
            return False

        checksum = 0
        for index2 in range(1,len(casrn)+1):
            checksum = checksum + int(casrn[-index2])*(index2-1)
        
        if str(checksum%10) != casrn[-1]:
            return False
        else:
            return True
 

    def sizeHint(self,option,index):
        """Generates sizehint for the view."""

        dispValue = index.data(QtCore.Qt.DisplayRole)
        chemical =  index.model().mapToSource(index).internalPointer().element
        #self.casrnXPath = etree.XPath('./CASNUMBER/DP_INTEGER1')
    
        if index.isValid():

            #get casrns.
            casrns = [x.text for x in self.casrnXPath(chemical) ]
            casrnstring = '<b>' + ', '.join(casrns) + '</b>'
            dispValue = casrnstring
    
            if index.isValid():
                align = QtCore.Qt.AlignLeft
                
            valueDoc = QtGui.QTextDocument()
            valueDoc.setTextWidth(option.rect.width())
            valueDoc.setHtml(dispValue)

            initialHeight = valueDoc.size().height()

            #valueDoc = self.addSvg(index,valueDoc)
            if index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgString')].toPyObject() != '':
                svgWidth = float(index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgWidth')].toPyObject())
                svgHeight = float(index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgHeight')].toPyObject())
                scaleFactor = valueDoc.textWidth()/svgWidth
                
                return  QtCore.QSize(option.rect.width(),initialHeight + svgHeight*scaleFactor+12)
            else:
                print index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgString')].toPyObject()
                return QtCore.QSize(option.rect.width(),initialHeight)

        else:
            return  QtCore.QSize(option.rect.width(),200)

    def addSvg(self,idx,doc):
        """Adds svg of chemical structure to the QtGui.QTextDocument. 

        This has been pre-made and is stored in the model."""

        #pprint(idx.data(QtCore.Qt.UserRole).toMap())
        svgString = idx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgString')].toPyObject()
        svgWidth = float(idx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgWidth')].toPyObject())
        svgHeight = float(idx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgHeight')].toPyObject())


        #note: for speed, the model contains the intrinsic width and height of 
        #the svg, as well as the svg. Unfortunately, the 
        #expensive part is the actual rendering, which must be done here.
        #This is done to speed up the real-time CAS-RN filtering. However,
        #I found the best way of doing this was actually to only start 
        #filtering when the length of the filter string is >2.
        #svgDict[smiles] has form [width,height,QSvgRenderer]
        scaleFactor = doc.textWidth()/svgWidth
        svgSize = QtCore.QSize(svgWidth,svgHeight)

        #make svg renderer
        myXmlStreamReader = QtCore.QXmlStreamReader(svgString)
        myRenderer = QtSvg.QSvgRenderer(myXmlStreamReader)


        myRenderer.setViewBox(QtCore.QRect(0,0,svgSize.width()*scaleFactor,svgSize.height()*scaleFactor))
        
        myImage = QtGui.QImage(svgSize*scaleFactor,QtGui.QImage.Format_ARGB32)
        
        myTempPainter = QtGui.QPainter(myImage)
        myTempPainter.fillRect(myImage.rect(),QtCore.Qt.white)
        #self.svgDict[smiles][2].render(myTempPainter)
        myRenderer.render(myTempPainter)
        myTempPainter.end()
        cursor = QtGui.QTextCursor(doc)
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText('\n')
        cursor.insertImage(myImage)
        return doc

class chemSFProxyModel(QtGui.QSortFilterProxyModel):
    """Lies between the model and the chemical view to enable filtering.
    
    This also makes sure only xml elements at the 'chemical' level are
    displayed."""

    def __init__(self):

        super(chemSFProxyModel,self).__init__()

        self.casrnFilterText = ''
        cs = QtCore.Qt.CaseInsensitive
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        self.blankRegExp = QtCore.QRegExp('',cs,syntax)
        self.setFilterRegExp(QtCore.QRegExp('',cs,syntax))

        #make a blank model for the boolean evaluator
        mdl = filterCollection()
        self.be = booleanEvaluator('',mdl)
        self.invalidateFilter()
        self.casrnXPath = etree.XPath('./CASNUMBER[ACTIVE[last()]/Q_QUAL_VALUE1=\'Y\']/DP_INTEGER1')

        self.sortMode = ''
        self.setDynamicSortFilter(True)
        self.sort(0)

    def filterAcceptsRow(self,row,parent):
        """Returns True if the row is to be visible, False otherwise."""
        #print self.rowCount(self.mapFromSource(parent))
    #if self.publishedCheckBox.isChecked():
        
        #print row,parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('xmlelement')].toPyObject(),parent.child(row,0).data(QtCore.Qt.UserRole).toMap()
        #print '---'
        #if the xml tag is 'CHEMICAL'
        if parent != QtCore.QModelIndex():
            if parent.child(row,0).internalPointer().element.tag == 'CHEMICAL':

                #print row,parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('xmlelement')].toPyObject().tag
                #get all CAS-RN's.
                chemXML = parent.child(row,0).internalPointer().element
                casrnElements = self.casrnXPath(chemXML)
                casrnElements = [x.text for x in casrnElements]
                #print casrnElements

                #casrnFilterResult = True
                #for casrn in casrnElements:
                #    
                #    if self.filterRegExp() != self.blankRegExp \
                #            and not self.filterRegExp().indexIn(casrn) >= 0:
                #        casrnFilterResult = False
                casrnFilterResult = False

                if self.filterRegExp() == self.blankRegExp:
                    casrnFilterResult = True
                else:
                    for casrn in casrnElements:
                        if self.filterRegExp().indexIn(casrn) >= 0:
                            casrnFilterResult = True
                            break
                            


                #print self.filterRegExp().indexIn(parent.child(row,0).data(QtCore.Qt.DisplayRole).toString()) >= 0
                #print casrnFilterResult

                #Now get the DPFilterWidget filter result.
                #DPFilterWidgetResult = True

                #get source index
                sourceIdx = parent.model().index(row,0,parent)

                DPFilterWidgetResult = self.be.applyToChemical(sourceIdx)

                return casrnFilterResult and DPFilterWidgetResult
            else:
                #not at the 'CHEMICAL' level
                return False
        else:
            return True

    def lessThan(self,leftIdx,rightIdx):
        """Order by CAS-RN."""

        leftCASRN = 0
        rightCASRN = 0
        if leftIdx != QtCore.QModelIndex():
            #get CASRNs
            leftChemXML = leftIdx.internalPointer().element
            leftCASRNElements = self.casrnXPath(leftChemXML)
            leftCASRNElements = [int(x.text) for x in leftCASRNElements]
        if leftCASRNElements != []:
            leftCASRN = min(leftCASRNElements)


        if rightIdx != QtCore.QModelIndex():
            rightChemXML = rightIdx.internalPointer().element
            rightCASRNElements = self.casrnXPath(rightChemXML)
            rightCASRNElements = [int(x.text) for x in rightCASRNElements]
        if rightCASRNElements != []:
            rightCASRN = min(rightCASRNElements)

        if leftCASRN < rightCASRN:
            return True
        else:
            return False
    def setAdvancedFilter(self,myBooleanEvaluator):
        """Applies the filter."""
        self.be = myBooleanEvaluator
        self.invalidateFilter()



class dataPointSFProxyModel(QtGui.QSortFilterProxyModel):
    """Lies between the model and the datapoint view to enable filtering."""
    
    def __init__(self,structureDict):

        super(dataPointSFProxyModel,self).__init__()

        self.structureDict = structureDict
        self.dataPointFilterBehaviour = 'Name'
        self.dataPointSortBy = 'Name'

        self.dataPointMostRecentOnly = False
        #self.setFilterKeyColumn(0)
        

        cs = QtCore.Qt.CaseInsensitive
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        self.setFilterRegExp(QtCore.QRegExp('',cs,syntax))


        self.reportOrder = ['CAS-RN',
            'Common name',
            'AICS Name',
            'RQ',
            'Release Mitigation Factor',
            'Volume (T) (8 May 2013)',
            'Mitigated PEC(river) (ug/L)',
            'PNEC (ug/L)',
            'Persistent?',
            'Bioaccumulative?',
            'Toxic?',
            'Comparison with Canada DSL',
            'Assessment conclusion',
            'Comments',
            'Internal Notes (not for publication)',
            'Group',
            'Assessment Status',
            'Batch no.',
            'Assessor',
            'Reasonable worst case exposure scenario appropriate?',
            'AHVICL Uses (2006)',
            'AHVICL Threshold Range (2006)',
            'SPIN UC62 - Total Tonnage',
            'SPIN Uses (UC62)',
            'Detailed SPIN Uses',
            'NICNAS ID',
            'IMAP Status',
            'Chemical class',
            'SMILES String',
            'Molecular formula',
            'Molecular weight (g/mol)',
            'High concern use for environment?',
            'Perfluorinated? (NICNAS master list 12 Oct 2012)',
            'Montreal',
            'SGG',
            'Rotterdam',
            'Stockholm',
            'REACH (SVHCs)',
            'EDC (US EPA)',
            'EDC (Europe)',
            'On DSL',
            'DSL P',
            'DSL B',
            'DSL iT',
            'log Kow',
            'Water solubility  (mg/L)',
            'Melting point (deg C)',
            'Vapour pressure (Pa)',
            'pKa',
            'Ionisable in the environment?',
            'Phys Chem Notes',
            'Obs. BOD (301C)',
            'BOD (301C)',
            'Primary halflife in water (days; Catalogic 301C)',
            'Reasons for Categorisation (P)',
            'Reasons for Categorisation (B)',
            'Toxicity Notes',
            'Fish ECOSAR (mg/L; Neutral Organics SAR)',
            'Daphnia ECOSAR (mg/L; Neutral Organics SAR)',
            'Algae ECOSAR (mg/L; Neutral Organics SAR)',
            'ECOSAR Acute fish endpoint (mg/L)',
            'ECOSAR Acute invertebrate endpoint (mg/L)',
            'ECOSAR Acute algae endpoint (mg/L)',
            'Fish endpoint (mg/L)',
            'Invertebrate endpoint (mg/L)',
            'Algae endpoint (mg/L)',
            'Other endpoint (mg/L)',
            'Pivotal endpoint (mg/L)',
            'Pivotal endpoint type',
            'Assessment factor',
            'AF Notes',
            'Default volume? (8 May 2013)',
            'PEC(river) (ug/L ; 8 May 2013)']

        self.invalidateFilter()
        self.invalidate()
        self.setDynamicSortFilter(True)
        self.sort(0)


    def filterAcceptsRow(self,row,parent):
        """
        Returns True if the row is to be visible, False otherwise.
        """
        #print self.rowCount(self.mapFromSource(parent))
        #print parent==QtCore.QModelIndex()
        #print self.dataPointMostRecentOnly
        #print row,parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('xmlelement')].toPyObject(),parent.child(row,0).data(QtCore.Qt.UserRole).toMap()
        #print '---'
        
        if parent != QtCore.QModelIndex():
            #if it's a DATA_POINT
            if QtCore.QString('DP_L1_KINGDOM_CODE') in parent.child(row,0)\
                .data(QtCore.Qt.UserRole).toMap().keys():
                #set all filters to True.
                mostRecent = True
                filterByName = True

                #modify filters according to what has been selected
                if self.dataPointMostRecentOnly == True:
                    if str(parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('CREATED_DATE')].toString()) != 'AUTO':
                            mostRecent = False

                if self.dataPointFilterBehaviour == 'Name':
                    #dFilter = (self.filterRegExp().indexIn(parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('DP_L0_DOMAIN_CODE')].toString()) >= 0 )
                    dFilter = (self.filterRegExp().indexIn(self.structureDict[str(parent.child(row,0).data(QtCore.Qt.UserRole).toMap()[QtCore.QString('DP_L0_DOMAIN_CODE')].toString())]['longName']) >= 0 )
                elif self.dataPointFilterBehaviour == 'Name and content':
                    dFilter = (self.filterRegExp().indexIn(parent.child(row,0).data(QtCore.Qt.DisplayRole).toString()) >= 0 )
               
                return (mostRecent and dFilter)
            
            else: # if it's not a DATA_POINT
                return True
            
        else:
            return True

    def lessThan(self,leftIdx,rightIdx):
        """Implements ordering of data points."""
        #print 'lessThan:',self.dataPointSortBy
        if self.dataPointSortBy == 'Name':
            if leftIdx.data(QtCore.Qt.DisplayRole).toString() < rightIdx.data(QtCore.Qt.DisplayRole).toString():
                return True
            else:
                return False

        elif self.dataPointSortBy == 'Report Order':

            leftXML = leftIdx.internalPointer().element
            rightXML = rightIdx.internalPointer().element
            if leftXML.tag not in self.structureDict or rightXML.tag not in self.structureDict:
                return True
            leftName = self.structureDict[leftXML.tag]['longName']
            rightName = self.structureDict[rightXML.tag]['longName']

            if leftName in self.reportOrder and rightName in self.reportOrder:
                leftOrderNum = self.reportOrder.index(leftName)
                rightOrderNum = self.reportOrder.index(rightName)
                #print [leftName], [rightName],leftOrderNum < rightOrderNum

                return leftOrderNum < rightOrderNum
            elif leftName in self.reportOrder and rightName not in self.reportOrder:
                return True
            elif leftName not in self.reportOrder and rightName in self.reportOrder:
                return False
            else:
                return True
        elif self.dataPointSortBy == 'Date':
            #xml is in date order.
            return True

        else:
            print 'lessThan in dataPointSFProxyModel: self.dataPointSortBy not one of accepted values.'
            print self.dataPointSortBy
            exit()



class qualifierSFProxyModel(QtGui.QSortFilterProxyModel):
    """Lies between the model and the qualifier view to enable filtering. 
    
    Currently does no sorting or filtering (Identity transformation)"""

    def __init__(self):

        super(qualifierSFProxyModel,self).__init__()
        
    def filterAcceptsRow(self,row,parent):
        """Returns True if the row is to be visible, False otherwise."""

        if parent != QtCore.QModelIndex():
            #if it's a QUALIFIER
            if QtCore.QString('Q_L1_KINGDOM_CODE') in \
                parent.child(row,0).data(QtCore.Qt.UserRole).toMap().keys():
                return True
            elif QtCore.QString('DP_L1_KINGDOM_CODE') in \
                parent.child(row,0).data(QtCore.Qt.UserRole).toMap().keys():
                #It appears to be necessary to be have the filter allow the parent.
                return True
            else:
                return False
        else:
            return False




class TreeItem(object):
    """Objects used by the underlying model. 
    
    Each TreeItem corresponds to an xml element in the underlying xml model.
    This is not a subclass of a PyQt class, but is used by the #TreeModel class
    (which is a subclass of QtCore.QAbstractItemModel). This allows access
    to sub-elements of the xml, using xpath. This is necessary to do things
    like view chemicals with a particular CAS-RN. "Chemical" TreeItems can be
    displayed if the underlying xml element has a subelement 
    ./CASNUMBER/DP_INTEGER1 which matches the CAS-RN in question, for example."""

    def __init__(self, element, parent=None):
        #print 'TreeItem: ',element,parent
        
        #save xml element for later access
        self.element = element

        #parent item is a TreeItem object
        self.parentItem = weakref.ref(parent) if parent else None

        #generate data dict from xml element
        leaves = element.xpath('./*[not(*)]')

        self.datadict = {}
        for item in leaves:
            self.datadict[item.tag] = item.text
    
        #try pass xml element to delegate. This returns a QVariant::UserType
        #object, which can be turned back into the etree element with 
        # toPyObject().
        #self.datadict['xmlelement']=self.element
        self.datadict['DP_L0_DOMAIN_CODE']=str(self.element.tag)
        #print str(self.element.tag)


        #self.datadict = QtCore.QVariant(self.datadict).toMap()
        self.childItems = []

    def appendChild(self, item):
        """Appends item to self's list of children."""

        self.childItems.append(item)

    def child(self, row):
        """Returns child number "row" of self. """
        try:
            return self.childItems[row]
        except:
            return None

    def childCount(self):
        """Returns number of children of self."""
        return len(self.childItems)

    def columnCount(self):
        """Returns 1. Only one column in this model."""
        return 1

    def data(self, column):
        """Returns datadict. 
        
        datadict contains the xml element corresponding to self 
        (key "xmlelement"), as well as the data stored in this xml
        element (keys DP_L0_DOMAIN_CODE, DP_QUAL_VALUE1 etc.)"""
        try:
            #return self.itemData[column]
            return self.datadict
        except IndexError:
            return None

    def parent(self):
        """Returns the parent TreeItem."""
        return (self.parentItem() if self.parentItem else None)

    def row(self):
        """Returns the index of the TreeItem in its parent's child list."""
        if self.parentItem:
            #note that this index function acts on a list.
            return self.parentItem().childItems.index(self)

        return 0

    def getSvg(self):
        """Returns the SVG for the chemical, if it exists."""
        try:
            svg = self.datadict[QtCore.QString('svgString')]
            if svg != '':
                return svg
            else:
                return None
        except:
            return None

    def getXmlElement(self):
        """Returns the xml element for the TreeItem."""
        return self.element




class TreeModel(QtCore.QAbstractItemModel):
    """The model used by the PyQt views.
    
    This acts as an interface to the xml version of the database. The 
    model is a tree of TreeItem objects, each of which has a list of 
    child TreeItems and a parent TreeItem. Three proxy models interface to
    this model and tell the three listviews (chemView, dataPointView and 
    qualifierView) what to display.
    """
    def __init__(self, rootelement,structureDict,parent=None):
        super(TreeModel, self).__init__(parent)

        #self.rootItem = TreeItem(etree.Element('aaa'))
    #self.rootItem.element.append(rootelement)
        ####self.rootItem = TreeItem(rootelement)
    ####
    #####print 'TreeModel init: ',rootelement,parent
        ####telement=self.treestep(rootelement,self.rootItem)
        #####telement=self.treestep(rootelement,parent)
        self.structureDict = structureDict
        print 'Making TreeModel'
        self.rootItem = self.treestep(rootelement,None)
        print 'done'
        #print 'self.rootItem tag ',self.rootItem.element[0].tag
        #print etree.tostring(rootelement,pretty_print=True)
        print 'running self.makeSvgs'

        self.makeSvgs()
        print 'done'

        
    

    def columnCount(self, parent):
        """Returns number of columns. (Always 1)"""
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        """Returns data depending on value of role.
        
        If role == QtCore.Qt.DisplayRole, this function returns a string
        for the view to display.

        If role == QtCore.Qt.UserRole, this function returns the datadict of
        the TreeItem corresponding to index."""


        dispValue = ''
        if not index.isValid():
            return None

        if role == QtCore.Qt.UserRole:
            item = index.internalPointer()
            return item.data(index.column())
        elif role == QtCore.Qt.DisplayRole:
            #get dict out of item.data for role QtCore.Qt.UserRole
            value = self.data(index,QtCore.Qt.UserRole)

            #print 'paint:',{str(key) : str(value.toMap()[key]) for key in value.toMap().keys()}
   
            dispValue = '<b>' +\
                value['DP_L0_DOMAIN_CODE'] + '</b> '
            
            if 'DP_L1_KINGDOM_CODE' in value.keys():
                #data point:
                dispValue = '<b>' +\
                    self.structureDict[value['DP_L0_DOMAIN_CODE']]\
                    ['longName'] +'</b> '

                DP_L1_KINGDOM_CODE = value['DP_L1_KINGDOM_CODE']

                if DP_L1_KINGDOM_CODE == 'DP_YN' :
                    dispValue = dispValue +  \
                            value['DP_QUAL_VALUE1']

                elif DP_L1_KINGDOM_CODE == 'DP_CLOB':
                    DP_CLOB1 = value['DP_CLOB1']

                    #basic sanitization
                    DP_CLOB1 = DP_CLOB1.replace('<','&lt;')
                    DP_CLOB1 = DP_CLOB1.replace('>','&gt;')

                    dispValue = dispValue + DP_CLOB1

                elif DP_L1_KINGDOM_CODE == 'DP_VARCHAR2' or \
                        DP_L1_KINGDOM_CODE == 'DP_VARCHAR2CHOICE':
                    DP_VARCHAR2 = value['DP_QUAL_VALUE1']

                    #basic sanitization
                    DP_VARCHAR2 = DP_VARCHAR2.replace('<','&lt;')
                    DP_VARCHAR2 = DP_VARCHAR2.replace('>','&gt;')

                    dispValue = dispValue + DP_VARCHAR2

                elif DP_L1_KINGDOM_CODE == 'DP_INTEGER':
                    dispValue = dispValue +  \
                            value['DP_INTEGER1']

                elif DP_L1_KINGDOM_CODE == 'DP_IMAPSTATUS':
                    dispValue = dispValue +  \
                            value['DP_QUAL_VALUE1']

                elif DP_L1_KINGDOM_CODE == 'DP_FLOAT':

                    #work out what sort of float we have.
                    if 'DP_QUANT_VALUE1' in value.keys() \
                        and  'DP_QUANT_VALUE2' not in value.keys() \
                        and  'DP_QUAL_VALUE1' not in value.keys() \
                        and  'DP_QUAL_VALUE2' not in value.keys():

                        #single float.
                        dispValue = dispValue +  \
                                value['DP_QUANT_VALUE1']
                        if 'DP_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['DP_UNIT_CODE1']
                            
                    elif 'DP_QUANT_VALUE1' in value.keys() \
                        and  'DP_QUAL_VALUE1' in value.keys() \
                        and  'DP_QUANT_VALUE2' not in value.keys() \
                        and  'DP_QUAL_VALUE2' not in value.keys():
                        
                        # greater than or greater than or equal to
                        DP_QUAL_VALUE1 =  value['DP_QUAL_VALUE1']
                        DP_QUAL_VALUE1=DP_QUAL_VALUE1.replace('>','&gt;')
                        dispValue = dispValue +  DP_QUAL_VALUE1 + \
                                value['DP_QUANT_VALUE1']
                        if 'DP_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['DP_UNIT_CODE1']
                    
                    elif 'DP_QUANT_VALUE2' in value.keys() \
                        and  'DP_QUAL_VALUE2' in value.keys() \
                        and  'DP_QUANT_VALUE1' not in value.keys() \
                        and  'DP_QUAL_VALUE1' not in value.keys():
                        # less than or less than or equal to
                        DP_QUAL_VALUE2 =  value['DP_QUAL_VALUE2']
                        DP_QUAL_VALUE2 = DP_QUAL_VALUE2.replace('<','&lt;')
                        dispValue = dispValue + DP_QUAL_VALUE2 +\
                                value['DP_QUANT_VALUE2']
                        if 'DP_UNIT_CODE2' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['DP_UNIT_CODE2']
                        
                    
                    elif 'DP_QUANT_VALUE1' in value.keys() \
                        and  'DP_QUAL_VALUE1' in value.keys() \
                        and  'DP_QUANT_VALUE2' in value.keys() \
                        and  'DP_QUAL_VALUE2' in value.keys():
                    
                        # range
                        if value['DP_QUAL_VALUE1'] == '>':
                            l = '('
                        elif value['DP_QUAL_VALUE1']=='>=':
                            l = '['
                        else:
                            dispValue = dispValue +' Unknown type of float ('
                            return dispValue

                        if value['DP_QUAL_VALUE2'] == '<':
                            r = ')'
                        elif value['DP_QUAL_VALUE2']=='<=':
                            r = ']'
                        else:
                            dispValue = dispValue +' Unknown type of float )'
                            return dispValue
                        
                        dispValue = dispValue + l+ \
                            value['DP_QUANT_VALUE1']\
                            +', '+\
                            value['DP_QUANT_VALUE2']+r
                        if 'DP_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['DP_UNIT_CODE1']
                        
                        if 'DP_UNIT_CODE2' in value.keys() and\
                            value['DP_UNIT_CODE1'] != \
                            value['DP_UNIT_CODE2']:
                            #tack DP_UNIT_CODE2 onto the end of dispValue
                            #if DP_UNIT_CODE2 != DP_UNIT_CODE1
                            dispValue = dispValue + ' ' + \
                                value['DP_UNIT_CODE2']
                    elif 'DP_QUAL_VALUE1' in value.keys() \
                            and 'DP_QUANT_VALUE1' not in value.keys()\
                            and 'DP_QUANT_VALUE2' not in value.keys():
                                dispValue = dispValue + ' ' + \
                                value['DP_QUAL_VALUE1']

            elif 'Q_L1_KINGDOM_CODE' in value.keys():
                #Qualifier.

                #get parent data.
                parentValue = self.data(index.parent(),QtCore.Qt.UserRole)
                parentDc = parentValue['DP_L0_DOMAIN_CODE']
                #dispValue = dispValue + ' aaa.'
                dc = value['DP_L0_DOMAIN_CODE']
                dispValue = '<b>' +\
                    self.structureDict[parentDc]['qualifiers'][dc]['longName'] +'</b> '

                Q_L1_KINGDOM_CODE = value['Q_L1_KINGDOM_CODE']

                if Q_L1_KINGDOM_CODE == 'Q_YN' or Q_L1_KINGDOM_CODE == 'Q_YNPStar':
                    dispValue = dispValue +  \
                            value['Q_QUAL_VALUE1']

                elif Q_L1_KINGDOM_CODE == 'Q_CLOB':
                    Q_CLOB1 = value['Q_CLOB1']

                    #basic sanitization
                    Q_CLOB1 = Q_CLOB1.replace('<','&lt;')
                    Q_CLOB1 = Q_CLOB1.replace('>','&gt;')

                    dispValue = dispValue + Q_CLOB1

                elif Q_L1_KINGDOM_CODE == 'Q_VARCHAR2':
                    Q_VARCHAR2 = value['Q_QUAL_VALUE1']

                    #basic sanitization
                    Q_VARCHAR2 = Q_VARCHAR2.replace('<','&lt;')
                    Q_VARCHAR2 = Q_VARCHAR2.replace('>','&gt;')

                    dispValue = dispValue + Q_VARCHAR2

                elif Q_L1_KINGDOM_CODE == 'Q_INTEGER':
                    dispValue = dispValue +  \
                            value['Q_INTEGER1']

                elif Q_L1_KINGDOM_CODE == 'Q_FLOAT':

                    #work out what sort of float we have.
                    if 'Q_QUANT_VALUE1' in value.keys() \
                        and  'Q_QUANT_VALUE2' not in value.keys() \
                        and  'Q_QUAL_VALUE1' not in value.keys() \
                        and  'Q_QUAL_VALUE2' not in value.keys():

                        #single float.
                        dispValue = dispValue +  \
                                value['Q_QUANT_VALUE1']
                        if 'Q_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['Q_UNIT_CODE1']
                            
                    elif 'Q_QUANT_VALUE1' in value.keys() \
                        and  'Q_QUAL_VALUE1' in value.keys() \
                        and  'Q_QUANT_VALUE2' not in value.keys() \
                        and  'Q_QUAL_VALUE2' not in value.keys():
                        
                        # greater than or greater than or equal to
                        Q_QUAL_VALUE1 =  value['Q_QUAL_VALUE1']
                        Q_QUAL_VALUE1.replace('>','&gt;')
                        dispValue = dispValue +  Q_QUAL_VALUE1 + \
                                value['Q_QUANT_VALUE1']
                        if 'Q_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['Q_UNIT_CODE1']
                    
                    elif 'Q_QUANT_VALUE2' in value.keys() \
                        and  'Q_QUAL_VALUE2' in value.keys() \
                        and  'Q_QUANT_VALUE1' not in value.keys() \
                        and  'Q_QUAL_VALUE1' not in value.keys():
                        # less than or less than or equal to
                        Q_QUAL_VALUE2 =  value['Q_QUAL_VALUE2']
                        Q_QUAL_VALUE2.replace('<','&lt;')
                        dispValue = dispValue + Q_QUAL_VALUE2 +\
                                value['Q_QUANT_VALUE2']
                        if 'Q_UNIT_CODE2' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['Q_UNIT_CODE2']
                        
                    
                    elif 'Q_QUANT_VALUE1' in value.keys() \
                        and  'Q_QUAL_VALUE1' in value.keys() \
                        and  'Q_QUANT_VALUE2' in value.keys() \
                        and  'Q_QUAL_VALUE2' in value.keys():
                    
                        # range
                        if value['Q_QUAL_VALUE1'] == '>':
                            l = '('
                        elif value['Q_QUAL_VALUE1']=='>=':
                            l = '['
                        else:
                            dispValue = dispValue +' Unknown type of float ('
                            return dispValue

                        if value['Q_QUAL_VALUE2'] == '<':
                            r = ')'
                        elif value['Q_QUAL_VALUE2']=='<=':
                            r = ']'
                        else:
                            dispValue = dispValue +' Unknown type of float )'
                            return dispValue
                        
                        dispValue = dispValue + l+ \
                            value['Q_QUANT_VALUE1']\
                            +', '+\
                            value['Q_QUANT_VALUE2']+r
                        if 'Q_UNIT_CODE1' in value.keys():
                            #tack unit onto the end of dispValue
                            dispValue = dispValue + ' ' + \
                                value['Q_UNIT_CODE1']
                        
                        if 'Q_UNIT_CODE2' in value.keys() and\
                            value['Q_UNIT_CODE1'] != \
                            value['Q_UNIT_CODE2']:
                            #tack DP_UNIT_CODE2 onto the end of dispValue
                            #if DP_UNIT_CODE2 != DP_UNIT_CODE1
                            dispValue = dispValue + ' ' + \
                                value['Q_UNIT_CODE2']
            #cb = index.internalPointer().element.xpath('./CREATED_BY/text()')[0]
            #cd = index.internalPointer().element.xpath('./CREATED_DATE/text()')[0]

            #return dispValue + ' '+cb + ' ' + cd
            return dispValue

    
    def setData(self,index,inDict,role):
        """Sets data at index to data contained in inDict for the given role."""

        def insertXml(dp,name,inDict,insertPosition):
            if name in inDict.keys():
                a = etree.Element(name)
                a.text = inDict[name]
                dp.insert(insertPosition,a)
                insertPosition += 1
            return insertPosition
        print 'setData: '

        parentIdx = index.parent()

        # get treeItem of interest, and its parent treeItem.
        ti = index.internalPointer()
        tiParent = ti.parent()

        #save some stuff for later.
        if 'DP_L1_KINGDOM_CODE' in inDict.keys():
            #it's a datapoint
            parentChemicalId = str(parentIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('CHEMICAL_ID')].toPyObject())

            dataPointId = str(index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('DATA_POINT_ID')].toPyObject())
            
            #check: Don't edit if DATA_POINT_ID_ID is not AUTO
            if dataPointId != 'AUTO':
                return
        elif 'Q_L1_KINGDOM_CODE' in inDict.keys():
            parentDataPointId = str(parentIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('DATA_POINT_ID')].toPyObject())
            qualifierId = str(index.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('QUALIFIER_ID')].toPyObject())
            
            #check: Don't edit if QUALIFIER_ID is not AUTO
            if qualifierId != 'AUTO':
                return

        else:
            print 'setData: not qualifier or datapoint!'
            exit()
        
        del ti.datadict
        ti.datadict = {}


        #remove leaves in current xml element (ti.element)
        leafpath = etree.ETXPath('./*[not(*)]')
        leaves = leafpath(ti.element)

        for item  in leaves:
            item.getparent().remove(item)
        
        if 'DP_L1_KINGDOM_CODE' in inDict.keys():
            #it's a datapoint

            #Fix the current xml element based on what's in inDict.
            DATA_POINT_ID = etree.Element('DATA_POINT_ID')
            DATA_POINT_ID.text = 'AUTO'
            ti.element.insert(0,DATA_POINT_ID)

            CHEMICAL_ID = etree.Element('CHEMICAL_ID')
            CHEMICAL_ID.text = parentChemicalId
            ti.element.insert(1,CHEMICAL_ID)

            insertPosition = insertXml(ti.element,'DP_L1_KINGDOM_CODE',inDict,2)

            CREATED_BY = etree.Element('CREATED_BY')
            CREATED_BY.text = 'AUTO'
            ti.element.insert(insertPosition,CREATED_BY)
            insertPosition += 1

            CREATED_DATE = etree.Element('CREATED_DATE')
            CREATED_DATE.text = 'AUTO'
            ti.element.insert(insertPosition,CREATED_DATE)
            insertPosition += 1

            eList = ['DP_QUANT_VALUE1','DP_QUANT_VALUE2',\
                    'DP_INTEGER1', 'DP_INTEGER2',\
                    'DP_CLOB1',\
                    'DP_CONFIDENTIAL1_YN', 'DP_CONFIDENTIAL2_YN',\
                    'DP_QUAL_VALUE1','DP_QUAL_VALUE2',\
                    'DP_UNIT_CODE1','DP_UNIT_CODE2']
            for name in eList:
                print insertPosition,[x for x in ti.element]
                insertPosition = insertXml(ti.element,name,inDict,insertPosition)


            print 'setData: ', etree.tostring(ti.element,pretty_print=True)

            #now set up ti.datadict (from treeItem.__init__)
            leaves = leafpath(ti.element)

            ti.datadict = {}
            for item in leaves:
                ti.datadict[item.tag] = item.text

            #ti.datadict['xmlelement']=ti.element
            ti.datadict['DP_L0_DOMAIN_CODE']=str(ti.element.tag)
            #print str(self.element.tag)


            #ti.datadict = QtCore.QVariant(ti.datadict).toMap()


            print 'setData: ',ti.datadict


        elif 'DP_L0_DOMAIN_CODE' in inDict.keys():
            #it's a qualifier

            #Fix the current xml element based on what's in inDict.
            QUALIFIER_ID = etree.Element('QUALIFIER_ID')
            QUALIFIER_ID.text = 'AUTO'
            ti.element.insert(0,QUALIFIER_ID)

            DATA_POINT_ID = etree.Element('DATA_POINT_ID')
            DATA_POINT_ID.text = parentDataPointId
            ti.element.insert(1,DATA_POINT_ID)


            insertPosition = insertXml(ti.element,'Q_L1_KINGDOM_CODE',inDict,2)
            
            CREATED_BY = etree.Element('CREATED_BY')
            CREATED_BY.text = 'AUTO'
            ti.element.insert(insertPosition,CREATED_BY)
            insertPosition += 1

            CREATED_DATE = etree.Element('CREATED_DATE')
            CREATED_DATE.text = 'AUTO'
            ti.element.insert(insertPosition,CREATED_DATE)
            insertPosition += 1

            eList = ['Q_QUANT_VALUE1','Q_QUANT_VALUE2',\
                    'Q_INTEGER1', 'Q_INTEGER2',\
                    'Q_CLOB1','Q_CONFIDENTIAL1_YN', 'Q_CONFIDENTIAL2_YN',\
                    'Q_QUAL_VALUE1','Q_QUAL_VALUE2',\
                    'Q_UNIT_CODE1','Q_UNIT_CODE2']
            
            for name in eList:
                print insertPosition,[x for x in ti.element]
                insertPosition = insertXml(ti.element,name,inDict,insertPosition)


            print 'setData: ', etree.tostring(ti.element,pretty_print=True)

            #now set up ti.datadict (from treeItem.__init__)
            leaves = leafpath(ti.element)

            ti.datadict = {}
            for item in leaves:
                ti.datadict[item.tag] = item.text

            #ti.datadict['xmlelement']=ti.element
            ti.datadict['DP_L0_DOMAIN_CODE']=str(ti.element.tag)
            #print str(self.element.tag)


            #ti.datadict = QtCore.QVariant(ti.datadict).toMap()


            print 'setData: ',ti.datadict

        else:
            print 'setData: unknown data type.'
            exit()

        self.dataChanged.emit(index,index)



    def flags(self, index):
        """Returns that all items in the view are selectable and enabled."""
        if not index.isValid():
            return QtCore.Qt.NoItemFlags

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        """Returns header information. Not used currently."""
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.rootItem.data(section)

        return None

    def index(self, row, column, parent):
        """Returns index corresponding to TreeItem at (row,column) of parent."""
        if not self.hasIndex(row, column, parent):
        #print 'index: returning default.'
            return QtCore.QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        """Returns the index of the parent of the TreeItem at "index"."""
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        try:
            parentItem = childItem.parent()
        except:
            parentItem = None
        #if parentItem == None:
        #    print None,parentItem == self.rootItem
        #else:
        #    print parentItem.element.tag,parentItem == self.rootItem

        if parentItem == self.rootItem or parentItem == None:
            return QtCore.QModelIndex()
        #if parentItem == None:
        #    return QtCore.QModelIndex()
        #else:
        #    return self.createIndex(parentItem.row(), 0, parentItem)
        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent):
        """Returns number of rows of parent (where parent is a TreeItem)."""
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

    def treestep(self,telement,tparent=None):
        """Generates model from xml elementTree.
        
        This function recursively steps through all of the elements of
        telement and returns treeItems. When applied to the top element
        of the xml tree, the result is a set of TreeItems with parents
        and children corresponding to the parent and child xml elements.

        This function generates the structure which the model uses as the
        interface to the underlying xml database."""

        #print telement,tparent
        #if tparent != None:
        #    print telement.tag,tparent.element.tag
        #else:
        #    print telement.tag,None

        
        #make new item
        newitem = TreeItem(telement,tparent)
        if tparent != None:
            tparent.appendChild(newitem)
           
        #get elements which have children
        withchildpath = etree.ETXPath('./*[*]') 
        for tchild in withchildpath(telement):
            self.treestep(tchild,newitem)

        if tparent == None:
            return newitem


    def smiles2svg(self,inSmiles):
        """Generates an svg structure from a SMILES string."""
        #print 'smiles2svg'

        #Deal with smiles made up of multiple distinct structures
        smilesList = inSmiles.split('.')

        svgList = []
        obconv = openbabel.OBConversion()
        obconv.SetInAndOutFormats('smi','svg')
        for smiles in smilesList:
            mol = openbabel.OBMol()
            readOk=obconv.ReadString(mol,smiles)
            if readOk == False:
                svgList.append(self.errorSvg)
            else:
                svg = obconv.WriteString(mol)
                svgList.append(svg)
        
        if len(smilesList) != len(svgList):
            print 'smiles2svg: len(smilesList) != len(svgList)'
            exit()

        if len(svgList) == 1:
            return svgList[0]
        elif len(svgList) == 0:
            return None
        elif len(svgList) > 1:
            #print 'packing'
            #pack, with a rectangle around the constituent components.
            p = organiseSVG.packer(svgList)
            #print 'packed'
            outSvg = p.makeSVG('best',True)
            outSvg =  etree.tostring(outSvg,pretty_print=True)

            #outSvg = '<?xml version="1.0"?>\n'+outSvg
            return outSvg
        else:
            print 'smiles2svg: ???'
            exit()


    def makeSvgs(self):
        """Generates svg graphics for all chemicals in the model."""

        #make error svg
        self.errorSvg = """<?xml version="1.0"?>
<svg version="1.1" id="topsvg"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
xmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="200px" height="200px" viewBox="0 0 100 100">
<title>OBDepict</title>
<rect x="0" y="0" width="100" height="100" fill="white"/>
<text text-anchor="middle" font-size="6" fill ="black" font-family="sans-serif"
x="50" y="98" ></text>
<g transform="translate(0,0)">
<svg width="100" height="100" x="0" y="0" viewBox="0 0 80 80"
font-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="2"  stroke-linecap="round">
<text x="5" y="75" fill="red"  stroke="red" stroke-width="1" font-size="96" >X</text>
</svg>
</g>
</svg>"""

        if os.path.isfile('cachedStructures.pickled'):
            print 'cached Structures!'
            with open('cachedStructures.pickled','r') as pickleFile:
                self.svgDict = cPickle.load(pickleFile)
        else:
            self.svgDict = {}

        #pre-generate svg graphics for use in chemical delegate.
        
        noOfChemicals = self.rootItem.childCount()
        if noOfChemicals == 1:
            realRoot = self.rootItem.child(0)
            noOfChemicals = realRoot.childCount()
        else:
            print 'makeSvgs: no Chemicals!'
            print 'rootItem.getElement',self.rootItem.getXmlElement()
            print 'noOfChemicals: ',noOfChemicals
            print 'self.rootItem.childCount()',self.rootItem.childCount()
            print etree.tostring(self.rootItem.element,pretty_print=True)
            return

        
        for i in range(noOfChemicals):

            #print '#',i
            chemicalTi = realRoot.child(i)
            CASRNS = ', '.join(chemicalTi.element.xpath('./CASNUMBER/DP_INTEGER1/text()'))
            print '#',i,CASRNS

            
            self.makeOneSvg(chemicalTi)

        #pickle self.svgDict for later use
        with open('cachedStructures.pickled','w') as pickleFile:
            cPickle.dump(self.svgDict,pickleFile)

        #index of first chemical
        firstIdx = self.createIndex(0,0,realRoot)

        #index of last chemical
        lastIdx= self.createIndex(0,noOfChemicals-1,realRoot)

        self.dataChanged.emit(firstIdx,lastIdx)

    def makeOneSvg(self,chemicalTi):
        """Generates the svg structure for a TreeItem."""
        chemical = chemicalTi.element
        smilesXPath = etree.XPath('./SMILES[ACTIVE[last()]/Q_QUAL_VALUE1=\'Y\']/DP_CLOB1')
        smilesList = smilesXPath(chemical)

        #print [x.text for x in smilesList]

        if len(smilesList) == 1:
            smiles = smilesList[0]
            smilesText = smiles.text.strip()

            #Get or make the smiles.
            if smilesText in self.svgDict.keys():
                outSvg = self.svgDict[smilesText]
            else:
                outSvg = self.smiles2svg(smilesText)
            self.svgDict[smilesText] = outSvg

            if outSvg != None:
                #get real size
                svgMatch = re.match('.*viewBox.*viewBox=\"[ \t]*([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]*\".*$',outSvg,re.DOTALL)
                if svgMatch:
                    svgWidth = svgMatch.group(3)
                    svgHeight = svgMatch.group(4)
                else:
                    print 'makeSvgs: no match for width/height.'
                    exit()
            else:
                outSvg = ''
                svgWidth = ''
                svgHeight = ''

        elif len(smilesList) == 0:
            outSvg = ''
            svgWidth = ''
            svgHeight = ''
        else:
            svgList=[]
            for smiles in smilesList:
                smilesText = smiles.text.strip()
                if smilesText in self.svgDict.keys():
                    svgList.append(self.svgDict[smilesText])
                else:
                    svgList.append(self.smiles2svg(smilesText))
                    self.svgDict[smilesText] = svgList[-1]
            #svgList = [self.smiles2svg(smiles.text.strip()) for smiles in smilesList]
            #pack the svgs as well as possible
            p = organiseSVG.packer(svgList)
            outSvg = p.makeSVG('best',False)
            outSvg = etree.tostring(outSvg)
            svgWidth = p.bestWidth
            svgHeight = p.bestHeight

        #add information to datadict for the relevant treeitem.
        chemicalTi.datadict[QtCore.QString(u'svgString')]=outSvg
        chemicalTi.datadict[QtCore.QString(u'svgWidth')]=svgWidth
        chemicalTi.datadict[QtCore.QString(u'svgHeight')]=svgHeight

        outTuple = (outSvg,svgWidth,svgHeight)

        #if emit == True:
        #    #get index of chemical in chemical.parent.childItems list
        #    listIndex = chemicalTi.parent().childItems.index(chemicalTi)

        #    changedIdx= self.createIndex(0,listIndex,chemicalTi.parent())

        #    self.dataChanged.emit(changedIdx,changedIdx)

        return outTuple


class cascadingViewer(QtGui.QWidget):
    """The cascading viewer widget.
    
    A new cascadingViewer widget is created whenever an xml file is opened, and
    inserted into the assessmentTool instance. When the file is closed, 
    the cascadingViewer widget is destroyed."""

    filterApplied = QtCore.pyqtSignal()
    stopAutoSave = QtCore.pyqtSignal()
    def __init__(self,xmlin,schemain,xmlFileName,at):
            
        #call constructor of parent class
        super(cascadingViewer,self).__init__(at)
        print 'cv\'s parent: ',self.parentWidget()

        self.xmlFileName = xmlFileName
        self.at = weakref.ref(at)

        #call method that makes the model
        self.initData(xmlin,schemain)

        #call method that generates the UI.
        self.initUI()

        self.initAutoSave()

    def initData(self,xmlin,schemain):
        """Initialises the data used by the widget."""
        
        #generate dictionary from schema.
        #self.schemaXML = etree.parse('newschema.rng').getroot()

        

        self.schemaXML = schemain
        self.schemaRNG = etree.RelaxNG(self.schemaXML)
        self.dictFromSchema()
        
        #xmlin = etree.parse('temp2.xml')

        #This is a hack to get the filter models to work correctly.
        #I can't get things to work without having an extra dummy
        #level above the ROOT level.
        #xmlroot = xmlin.getroot()
        xmlroot = xmlin
        xmlrootminusone = etree.Element('ROOTMINUSONE')
        xmlrootminusone.append(xmlroot)

        self.model = TreeModel(xmlrootminusone,self.structureDict)

        #set up proxy model for filtering and sorting.

        self.chemProxyModel = chemSFProxyModel()
        self.chemProxyModel.dumpObjectInfo()
        #modelTester2 = ModelTest(self.chemProxyModel,self)
        self.chemProxyModel.setSourceModel(self.model)
        
        self.dataPointProxyModel = dataPointSFProxyModel(self.structureDict)
        self.dataPointProxyModel.setSourceModel(self.model)

        self.qualifierProxyModel = qualifierSFProxyModel()
        self.qualifierProxyModel.setSourceModel(self.model)

        #regexp = QtCore.QRegExp('CASNUMBER',QtCore.Qt.CaseInsensitive,QtCore.QRegExp.FixedString)
        #self.proxyModel.setFilterRegExp(regexp)
        

    
    def initUI(self):
        """Initialises the user interface."""

        self.chemView = QtGui.QListView()
        self.chemView.setAcceptDrops(True)
        #self.chemView.setAlternatingRowColors(True)
        self.chemView.setResizeMode(QtGui.QListView.Adjust)
        #self.chemView.setModel(self.proxyModel)
        self.chemView.setModel(self.chemProxyModel)
        #self.chemViewDe = MyDelegate()
        self.chemViewDe = MyChemicalDelegate()
        self.chemView.setItemDelegate(self.chemViewDe)
        self.chemSelectionModel = self.chemView.selectionModel()
        self.connect(self.chemSelectionModel,QtCore.SIGNAL('currentChanged(QModelIndex,QModelIndex)'),self.onChemicalViewChanged)

        self.dataPointView = QtGui.QListView()
        self.dataPointView.setAlternatingRowColors(True)
        self.dataPointView.setResizeMode(QtGui.QListView.Adjust)
        self.dataPointView.setModel(self.dataPointProxyModel)
        self.dataPointViewDe = MyDataPointDelegate()
        self.dataPointView.setItemDelegate(self.dataPointViewDe)
        self.dataPointSelectionModel = self.dataPointView.selectionModel()
        self.connect(self.dataPointSelectionModel,QtCore.SIGNAL('currentChanged(QModelIndex,QModelIndex)'),self.onDataPointViewChanged)
        


        self.qualifierView = QtGui.QListView()
        self.qualifierView.setAlternatingRowColors(True)
        self.qualifierView.setResizeMode(QtGui.QListView.Adjust)
        self.qualifierView.setModel(self.qualifierProxyModel)
        self.qualifierViewDe = MyQualifierDelegate()
        self.qualifierView.setItemDelegate(self.qualifierViewDe)
        self.qualifierSelectionModel = self.qualifierView.selectionModel()

        #self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem)))
        #self.chemView.setRootIndex(QtCore.QModelIndex())
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))

        #initialise selection
        #self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem)),QtGui.QItemSelectionModel.SelectCurrent)
        #self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0).child(0))),QtGui.QItemSelectionModel.SelectCurrent)
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem)).child(0,0),QtGui.QItemSelectionModel.SelectCurrent)
        #horizontal splitter that contains the three vertical layouts.
        hsplitter = QtGui.QSplitter()
        #self.onSplitterMoved emits resize signals of delegates to 
        #make the pictures resize properly.
        hsplitter.splitterMoved.connect(self.onSplitterMoved)

        #Create widget wrappers for vertical layouts
        chemWidget = QtGui.QWidget()
        dataPointWidget = QtGui.QWidget()
        qualifierWidget = QtGui.QWidget()

        chemVbox = QtGui.QVBoxLayout()
        dataPointVbox =  QtGui.QVBoxLayout()
        qualifierVbox = QtGui.QVBoxLayout()

        chemWidget.setLayout(chemVbox)
        dataPointWidget.setLayout(dataPointVbox)
        qualifierWidget.setLayout(qualifierVbox)

        #buttons for adding new chemical/datapoint/qualifier
        chemAddNewBtn = QtGui.QPushButton('New Chemical')
        dataPointAddNewBtn = QtGui.QPushButton('New Data Point')
        qualifierAddNewBtn = QtGui.QPushButton('New Qualifier')


        #set up datapoint sorting/filtering widget
        dataPointSortFilterForm = QtGui.QFormLayout()
        dataPointSortDPFilterWidget = QtGui.QWidget()
        dataPointSortDPFilterWidget.setLayout(dataPointSortFilterForm)

        #combo boxes/line edits for datapoint sorting/filtering widget
        self.dataPointFilterPatternLineEdit = QtGui.QLineEdit()
        self.dataPointFilterPatternLineEdit.textChanged.connect(self.onDataPointFilterPatternChanged)
        self.dataPointFilterBehaviourCombo = QtGui.QComboBox()
        self.dataPointFilterBehaviourCombo.addItem('Name')
        self.dataPointFilterBehaviourCombo.addItem('Name and content')
        for i in [0,1]:
            self.dataPointFilterBehaviourCombo.setItemData(i,QtCore.Qt.green,QtCore.Qt.BackgroundRole)

        self.dataPointSortByCombo = QtGui.QComboBox()
        self.dataPointSortByCombo.addItem('Report Order')
        self.dataPointSortByCombo.addItem('Name')
        self.dataPointSortByCombo.addItem('Date')
        self.dataPointSortByCombo.currentIndexChanged.connect(self.onDataPointSort)
        self.onDataPointSort()

        #select first of the (possibly sorted) datapoints. Only need to 
        #worry about this when the program is first opened - onDataPointViewChanged
        #takes care of it afterwards.
        currentlySelectedChemIdx = self.chemSelectionModel.currentIndex()
        currentlySelectedChemIdx = self.chemProxyModel.mapToSource(currentlySelectedChemIdx)
        currentlySelectedChemIdx = self.dataPointProxyModel.mapFromSource(currentlySelectedChemIdx)
        newDpIdx = self.dataPointProxyModel.index(0,0,currentlySelectedChemIdx)
        self.dataPointSelectionModel.setCurrentIndex(newDpIdx,QtGui.QItemSelectionModel.SelectCurrent)

        self.dataPointMostRecentOnlyCheckBox = QtGui.QCheckBox('Most recent only')
        self.dataPointMostRecentOnlyCheckBox.stateChanged.connect(self.onDataPointFilterPatternChanged)
        self.dataPointMostRecentOnlyCheckBox.setChecked(False)

        
        dataPointSortFilterForm.addRow('Filter Pattern:',self.dataPointFilterPatternLineEdit)
        dataPointSortFilterForm.addRow('Filter',self.dataPointFilterBehaviourCombo)
        dataPointSortFilterForm.addRow('Sort by:',self.dataPointSortByCombo)
        dataPointSortFilterForm.addRow(self.dataPointMostRecentOnlyCheckBox)

        #widget for sorting/filtering chemical view
        self.fw = DPFilterWidget(self.structureDict,self)
        print 'dp filter model',self.fw.mdl
        
        #Hook up the Set Filter button to the chemSFProxyModel. 
        #applyChemicalFilterToProxyModel is a wrapper that, strangely enough,
        #applies the chemical filter to the proxy model.
        self.fw.setFilterBtn.clicked.connect(self.applyChemicalFilterToProxyModel)
        self.fw.boolLineEdit.returnPressed.connect(self.applyChemicalFilterToProxyModel)

        self.showChemFiltersBtn = QtGui.QPushButton()
        #self.showChemFiltersBtn.setIcon(QtGui.QIcon.fromTheme('up',QtGui.QIcon('up.png')))
        self.showChemFiltersBtn.setIcon(self.style().standardIcon(QtGui.QStyle.SP_ArrowUp))
        self.hideChemFiltersBtn = QtGui.QPushButton()
        #self.hideChemFiltersBtn.setIcon(QtGui.QIcon.fromTheme('down',QtGui.QIcon('down.png')))
        self.hideChemFiltersBtn.setIcon(self.style().standardIcon(QtGui.QStyle.SP_ArrowDown))
        self.showChemFiltersBtn.clicked.connect(self.showChemFilters)
        self.hideChemFiltersBtn.clicked.connect(self.hideChemFilters)

        chemSortFilterForm = QtGui.QFormLayout()
        chemSortFilterVBox = QtGui.QVBoxLayout()
        chemSortFilterVBox.addWidget(self.showChemFiltersBtn)
        chemSortFilterVBox.addWidget(self.hideChemFiltersBtn)
        chemSortFilterVBox.addLayout(chemSortFilterForm)
        chemSortFilterVBox.addWidget(self.fw)
       
        self.chemSortDPFilterWidget = QtGui.QWidget()
        
        self.chemSortDPFilterWidget.setLayout(chemSortFilterVBox)

        self.CASRNFilterLineEdit = QtGui.QLineEdit()
        casrnValidator = intBlankValidator()
        self.CASRNFilterLineEdit.setValidator(casrnValidator)

        chemSortFilterForm.addRow('CAS-RN:',self.CASRNFilterLineEdit)
        
        #uncomment for live casrn filtering.
        #self.CASRNFilterLineEdit.textChanged.connect(self.onCASRNFilterPatternChanged)
        self.CASRNFilterLineEdit.returnPressed.connect(self.onCASRNFilterPatternChanged)
        #self.publishedCheckBox = QtGui.QCheckBox()
        #self.publishedCheckBox.setChecked(False)
        #chemSortFilterForm.addRow('Published',self.publishedCheckBox)

        #connect buttons
        chemAddNewBtn.clicked.connect(self.newChemical)
        dataPointAddNewBtn.clicked.connect(self.newDataPoint)
        qualifierAddNewBtn.clicked.connect(self.newQualifier)

        #set up editing
        self.dataPointView.doubleClicked.connect(self.editOrCopyDataPoint)
        self.qualifierView.doubleClicked.connect(self.editOrCopyQualifier)

        #set up context menus.
        self.makeContextMenus()

        #set up layouts
        chemVbox.addWidget(self.chemView)
        chemVbox.addWidget(self.showChemFiltersBtn)
        chemVbox.addWidget(self.hideChemFiltersBtn)
        chemVbox.addWidget(self.chemSortDPFilterWidget)
        chemVbox.addWidget(chemAddNewBtn)
        dataPointVbox.addWidget(self.dataPointView)
        dataPointVbox.addWidget(dataPointSortDPFilterWidget)
        dataPointVbox.addWidget(dataPointAddNewBtn)
        qualifierVbox.addWidget(self.qualifierView)
        qualifierVbox.addWidget(qualifierAddNewBtn)

        hsplitter.addWidget(chemWidget)
        hsplitter.addWidget(dataPointWidget)
        hsplitter.addWidget(qualifierWidget)
        aaa = QtGui.QVBoxLayout()
        aaa.addWidget(hsplitter)
        self.setLayout(aaa)
        
        self.hideChemFilters()
 
        self.setWindowTitle('Cascading Viewer')
        #self.setGeometry(50,50,800,700)


        self.show()

    def onSplitterMoved(self,pos,index):
        """Hack to make the list views resize correctly.
        
        This is called when moveSplitter is emitted from the splitter containing
        the three list views. It emits the sizeHintChanged signal of the
        delegates."""
        if index == 1:
            #splitter between chemical and data point views
            self.chemViewDe.sizeHintChanged.emit(QtCore.QModelIndex())
            self.dataPointViewDe.sizeHintChanged.emit(QtCore.QModelIndex())
        elif index == 2:
            #splitter between data point and qualifier views
            self.dataPointViewDe.sizeHintChanged.emit(QtCore.QModelIndex())
            self.qualifierViewDe.sizeHintChanged.emit(QtCore.QModelIndex())


    def onDataPointSort(self):
        """Re-sorts data points.
        
        Conneted to dataPointSortByCombo.currentIndexChanged."""
        print [self.dataPointSortByCombo.currentText()]
        self.dataPointProxyModel.dataPointSortBy = str(self.dataPointSortByCombo.currentText())
        self.dataPointProxyModel.invalidate()

        self.dataPointProxyModel.sort(0)



    def initAutoSave(self):
        """Initialises the auto-save thread.
        
        The autosaveWorker class that does all the work is found in autosave.py."""

        #set up autosave dir if it doesn't exist
        self.autoSavePath = os.path.join(os.getcwd(),'autosave')
        print self.autoSavePath
        if not os.path.isdir(self.autoSavePath):
            if os.path.exists(self.autoSavePath):
                print 'autosave: '+self.autoSavePath+'exists but isn\'t a directory'
                exit()
            else:
                #make it
                os.makedirs(self.autoSavePath)

        #Need to use old-style signals and slots here,
        #or disconnect slots manually. (Above new-style
        #signals and slots crashes on Windows.)
        self.autoSaveThread = QtCore.QThread(self)
        #self.autoSaveThread.setObjectName('autoSave')
        self.asw = autosave.autoSaveWorker(self,self.autoSavePath)
        self.asw.moveToThread(self.autoSaveThread)
        #TODO: connect error signal from asw up here.
        QtCore.QObject.connect(self.autoSaveThread,QtCore.SIGNAL('started()'),\
                self.asw.run)
        #stop autosave timer when stopAutoSave signal is fired
        QtCore.QObject.connect(self,QtCore.SIGNAL('stopAutoSave()'),\
                self.asw.killAutosaveWorker, QtCore.Qt.BlockingQueuedConnection)
        #quit the thread when autoSaveWorker.finished is fired
        #QtCore.QObject.connect(self.asw,QtCore.SIGNAL('finished()'),\
        #        self.autoSaveThread.quit)
        #mark the autosave worker for destruction when 
        #autoSaveWorker.finished is fired
        #QtCore.QObject.connect(self.asw,QtCore.SIGNAL('finished()'),\
        #        self.asw.deleteLater)
        #mark the thread for destruction when autoSaveThread.finished is fired.
        #QtCore.QObject.connect(self.autoSaveThread,QtCore.SIGNAL('finished()'),\
        #        self.autoSaveThread.deleteLater)
        #Connect up the slot so that when the autosave thread is finished,
        #cv is deleted.
        #QtCore.QObject.connect(self.autoSaveThread,QtCore.SIGNAL('destroyed()'),\
        #        self.at().cleanUpCv, QtCore.Qt.BlockingQueuedConnection)

        #hook up slot to do the saving
        QtCore.QObject.connect(self.asw,QtCore.SIGNAL('saveCurrentXml()'),\
                self.autoSaveXml)
        #save initial xml
        shutil.copy(self.xmlFileName,os.path.join(self.autoSavePath,'current.xml'))


        self.autoSaveThread.start()

        ##Switch off autosave for now, just connect up the signal to the
        ##cleanup function.
        #QtCore.QObject.connect(self,QtCore.SIGNAL('stopAutoSave()'),\
        #        self.at().cleanUpCv)


    def autoSaveXml(self):
        def tempsavexml(saveFileName):
                """Saves curent xml to saveFileName."""
            #try:
                #print saveFileName
                with open(saveFileName,'w') as saveFile:
                    sourceRoot = self.model.createIndex(0,0,self.model.rootItem).\
                            child(0,0).internalPointer().element
                    xml = copy.deepcopy(sourceRoot)

                    #print etree.tostring(tempXml,pretty_print=True)

                    #print len(xml)
                    #print xml

                    saveFile.write('<ROOT>\n')
                    for item in xml:
                        a = etree.tostring(item,pretty_print=True)
                        saveFile.write(a)
                    saveFile.write('</ROOT>\n')


            #except IOError:
            #        QtGui.QMessageBox.information(self,'Cascading Viewer',\
            #            'Couldn\'t save for some reason.',\
            #            QtGui.QMessageBox.Ok)


        maxIdx = 5

        maxPath = os.path.join(self.autoSavePath,str(maxIdx)+'.xml')
        if os.path.isfile(maxPath):
            os.remove(maxPath)

        for idx in reversed(range(1,maxIdx)):
            oldpath = os.path.join(self.autoSavePath,str(idx)+'.xml')
            newpath = os.path.join(self.autoSavePath,str(idx+1)+'.xml')
            
            if os.path.isfile(oldpath):
                shutil.move(oldpath,newpath)


        file1 = os.path.join(self.autoSavePath,'1.xml')
        tempsavexml(file1)


    def processAutoSaveError(self):
        """This will process auto-save errors when it is implemented."""
        errorStr = 'Error handling not implemented yet'
        print 'autosave: ',errorStr



    def showChemFilters(self):
        """Slot that makes the chemical filter interface visible."""
        self.showChemFiltersBtn.setVisible(False)
        self.hideChemFiltersBtn.setVisible(True)
        self.chemSortDPFilterWidget.setVisible(True)

    
    def hideChemFilters(self):
        """Slot that hides the chemical filter interface."""
        self.showChemFiltersBtn.setVisible(True)
        self.hideChemFiltersBtn.setVisible(False)
        self.chemSortDPFilterWidget.setVisible(False)
        
    def applyChemicalFilterToProxyModel(self):
        """Applies advanced filter to chemSFProxyModel.
        
        This slot is triggered when the \'Set filter\' button is pressed in the 
        DPFilterWidget. All of this is hooked up in initUI()."""
        #need this to make sure CAS-RN lineEdit is interrogated
        self.onCASRNFilterPatternChanged()

        checkString = self.fw.be.checkTokenizedString(self)
        if checkString  == True:
            self.chemProxyModel.setAdvancedFilter(self.fw.getBooleanEvaluator())
            self.filterApplied.emit()


    def onChemicalViewChanged(self,newIdx,oldIdx):
        """Resets the datapoint view."""
        
        newIdxBase = self.chemProxyModel.mapToSource(newIdx)

        newIdxNewModel = self.dataPointProxyModel.mapFromSource(newIdxBase)

        self.dataPointView.setRootIndex(newIdxNewModel)
        self.dataPointSelectionModel.setCurrentIndex(newIdxNewModel.child(0,0),QtGui.QItemSelectionModel.SelectCurrent)
    
    def onDataPointViewChanged(self,newIdx,oldIdx):
        """Resets the qualifier view."""
        
        newIdxBase = self.dataPointProxyModel.mapToSource(newIdx)

        newIdxNewModel = self.qualifierProxyModel.mapFromSource(newIdxBase)
        self.qualifierView.setRootIndex(newIdxNewModel)
        self.qualifierSelectionModel.setCurrentIndex(newIdxNewModel.child(0,0),QtGui.QItemSelectionModel.SelectCurrent)

    def onDataPointFilterPatternChanged(self):
        """Applies data point filtering."""

        #Reset the SFProxyModel
        self.dataPointProxyModel.dataPointFilterBehaviour = self.dataPointFilterBehaviourCombo.currentText()
        self.dataPointProxyModel.dataPointSortBy = self.dataPointSortByCombo.currentText()
        self.dataPointProxyModel.dataPointMostRecentOnly = self.dataPointMostRecentOnlyCheckBox.isChecked()

        #set regular expression
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        cs = QtCore.Qt.CaseInsensitive
        re = self.dataPointFilterPatternLineEdit.text()
        regExp = QtCore.QRegExp(re,cs,syntax)
        #just needed to trigger reapplication of filter.
        #self.proxyModel.filterRegExp is not used.
        self.dataPointProxyModel.setFilterRegExp(regExp)
        #print 'dp filter pattern changed'

    def onCASRNFilterPatternChanged(self):
        """Filters chemicals by CAS-RN."""
        print 'onCASRNFilterPatternChanged'
        syntax = QtCore.QRegExp.PatternSyntax(QtCore.QRegExp.RegExp)
        cs = QtCore.Qt.CaseInsensitive
        txt = str(self.CASRNFilterLineEdit.text())
        if txt == str(''):
            re = ''
        else:
            re = '^'+self.CASRNFilterLineEdit.text()
        regExp = QtCore.QRegExp(re,cs,syntax)
        print regExp

        self.chemProxyModel.setFilterRegExp(regExp)
        self.filterApplied.emit()
    
    def editOrCopyDataPoint(self,idx):
        """Edit data point or make a copy for editing. 
        
        If DATA_POINT_ID == "AUTO", really edit the datapoint, since it hasn't
        been uploaded to the remote database yet. Otherwise, make a copy of the 
        datapoint that points to the old one, so that there is always a record
        of the history of every piece of data."""
        #get DATA_POINT_ID. If it isn't 'AUTO', copy it. Otherwise,
        #edit it.

        print 'dpid: ', idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()
        dpid = idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'DATA_POINT_ID')].toString()
        print 'dpid: ',[dpid]

        if dpid == 'AUTO':
            self.editDataPoint(idx)
        else:
            #get datapoint name
            itemName = str(idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'DP_L0_DOMAIN_CODE')].toPyObject())
            chemIdx = idx.model().mapToSource(idx).parent()

            if self.dataPointAllowed(itemName,chemIdx):
                copyIdx = self.copyDataPoint(idx)
                self.editDataPoint(copyIdx)
    
    def editDataPoint(self,idx):
        """Edit data point. 
        
        Edit the datapoint, since it hasn't been uploaded to the 
        remote database yet"""

        print 'editDataPoint ',idx

        #get datapoint from underlying model
        dpIdx = idx.model().mapToSource(idx)
        print dpIdx
        chemti = dpIdx.internalPointer().parent()

        #make 'normal' dictionary (strings) from Qt types for input to
        #dataDialog.

        tempDict = {}
        for key in dpIdx.data(QtCore.Qt.UserRole).toMap().keys():
            tempDict[str(key)] = str(dpIdx.data(QtCore.Qt.UserRole).toMap()[key].toPyObject())


        #dg = dataDialog(self,underlyingIdx.data(QtCore.Qt.UserRole).toMap())
        dg = dataDialog(self,self.structureDict,tempDict)
        dgStatus = dg.exec_()
        if dgStatus:
            dgData = dg.getData()
            print '---'
            print dgData
            print '---'

            dpIdx.model().setData(dpIdx,dgData,QtCore.Qt.UserRole)
            if dgData['DP_L0_DOMAIN_CODE'] == 'SMILES':
                print 'Update figure!'
                self.model.makeOneSvg(chemti)
                #self.onSplitterMoved(1,QtCore.QModelIndex())
                self.chemViewDe.sizeHintChanged.emit(dpIdx.parent())


    def copyDataPoint(self,idx):
        """Copies a datapoint.
        
        Invoked when the user double clicks on a data point to edit it, and 
        the datapoint has been downloaded from the database. A copy of the 
        data point is created for the user to edit."""

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)
        print 'copyDataPoint: ',underlyingIdx

        ti = underlyingIdx.internalPointer()
        xmlCopy = copy.deepcopy(ti.element)
        
      
        chemidXPath = etree.XPath('.//CHEMICAL_ID')
        chemids = chemidXPath(xmlCopy)

        for item in chemids:
            item.text = 'AUTO'


        dpidXPath = etree.XPath('.//DATA_POINT_ID')
        dpids = dpidXPath(xmlCopy)

        for item in dpids:
            item.text = 'AUTO'


        qidXPath = etree.XPath('.//QUALIFIER_ID')
        qids = qidXPath(xmlCopy)

        for item in qids:
            item.text = 'AUTO'

        cbXPath = etree.XPath('.//CREATED_BY')
        cbs = cbXPath(xmlCopy)

        for item in cbs:
            item.text = 'AUTO'

        cdXPath = etree.XPath('.//CREATED_DATE')
        cds = cdXPath(xmlCopy)

        for item in cds:
            item.text = 'AUTO'
        print 'copyDataPoint: ',etree.tostring(xmlCopy,pretty_print=True)


        chemIdx = underlyingIdx.parent()
        dpIdx = underlyingIdx


        #add this to the appropriate element of the model. Need to
        #get the model index from the proxy model.
        
        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = chemIdx.internalPointer()
        dpti = dpIdx.internalPointer()

        #print 'chemti: ',chemti,chemti.parent()
        #print 'dpti: ',dpti

        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()
        chemti.element.append(xmlCopy)
        self.model.treestep(xmlCopy,chemti)
        #print etree.tostring(self.model.rootItem.element,pretty_print=True)
        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()

        #indexes have all changed. Reset the root indices of the views.
        #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
        #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
        #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
        #self.chemView.setRootIndex(QtCore.QModelIndex())
        #self.dataPointView.setRootIndex(QtCore.QModelIndex().child(0,0))
        #self.qualifierView.setRootIndex(QtCore.QModelIndex().child(0,0).child(0,0))
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))


        #set the current index of the chemical view to the 
        #index of chemti. Need to make new index in self.model and
        #map to self.proxyModel.
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
        #similarly, set current index of data point view
        #self.dataPointSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
        #set current index of data point view to newly created
        #data point.
        self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)

        print etree.tostring(xmlCopy,pretty_print=True)

        #get index of new datapoint

        outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1]))


        #dummy - want to return the dataPointProxyModel index for 
        #the new data point.

        return outIdx

    def copyDataPointToVisible(self,idx):
        """Copies a datapoint to all visible chemicals."""
        

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)

        ti = underlyingIdx.internalPointer()
        xmlCopy = copy.deepcopy(ti.element)
        
        chemidXPath = etree.XPath('.//CHEMICAL_ID')
        chemids = chemidXPath(xmlCopy)

        for item in chemids:
            item.text = 'AUTO'
        
        dpidXPath = etree.XPath('.//DATA_POINT_ID')
        dpids = dpidXPath(xmlCopy)

        for item in dpids:
            item.text = 'AUTO'


        qidXPath = etree.XPath('.//QUALIFIER_ID')
        qids = qidXPath(xmlCopy)

        for item in qids:
            item.text = 'AUTO'

        cbXPath = etree.XPath('.//CREATED_BY')
        cbs = cbXPath(xmlCopy)

        for item in cbs:
            item.text = 'AUTO'

        cdXPath = etree.XPath('.//CREATED_DATE')
        cds = cdXPath(xmlCopy)

        for item in cds:
            item.text = 'AUTO'
        print 'copyDataPoint: ',etree.tostring(xmlCopy,pretty_print=True)


        chemIdx = underlyingIdx.parent()
        dpIdx = underlyingIdx


        #add this to the appropriate element of the model. Need to
        #get the model index from the proxy model.
        
        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = chemIdx.internalPointer()
        dpti = dpIdx.internalPointer()

        #print 'chemti: ',chemti,chemti.parent()
        #print 'dpti: ',dpti

        sourceRoot = self.model.createIndex(0,0,self.model.rootItem).\
                child(0,0)
        proxyRoot = self.chemProxyModel.mapFromSource(sourceRoot)

        #for row in xrange(self.cv.chemProxyModel.rowCount(proxyRoot)):


        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()

        for row in xrange(self.chemProxyModel.rowCount(proxyRoot)):
            sourceVisTiIdx = proxyRoot.model().mapToSource(proxyRoot.child(row,0))
            visTi = sourceVisTiIdx.internalPointer()
            #copy to all chemicals except the selected chemical
            if visTi != chemti:
                xmlCopy2 = copy.deepcopy(xmlCopy)
                visTi.element.append(xmlCopy2)
                self.model.treestep(xmlCopy2,visTi)
        #print etree.tostring(self.model.rootItem.element,pretty_print=True)
        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()

        #indexes have all changed. Reset the root indices of the views.
        #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
        #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
        #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
        #self.chemView.setRootIndex(QtCore.QModelIndex())
        #self.dataPointView.setRootIndex(QtCore.QModelIndex().child(0,0))
        #self.qualifierView.setRootIndex(QtCore.QModelIndex().child(0,0).child(0,0))
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))


        #set the current index of the chemical view to the 
        #index of chemti. Need to make new index in self.model and
        #map to self.proxyModel.
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
        #similarly, set current index of data point view
        #self.dataPointSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
        #set current index of data point view to newly created
        #data point.
        self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)

        #print etree.tostring(xmlCopy,pretty_print=True)

        #get index of new datapoint

        #outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1]))

        outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti))


        #dummy - want to return the dataPointProxyModel index for 
        #the new data point.

        return outIdx

    def markAsInactive(self,idx):
        """Marks datapoint at idx as inactive."""
        

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)

        ti = underlyingIdx.internalPointer()
        xmlCopy = copy.deepcopy(ti.element)
        
        inactiveStr = """    <ACTIVE>
      <QUALIFIER_ID>AUTO</QUALIFIER_ID>
      <DATA_POINT_ID>AUTO</DATA_POINT_ID>
      <Q_L1_KINGDOM_CODE>Q_YN</Q_L1_KINGDOM_CODE>
      <CREATED_BY>AUTO</CREATED_BY>
      <CREATED_DATE>AUTO</CREATED_DATE>
      <Q_CONFIDENTIAL1_YN>N</Q_CONFIDENTIAL1_YN>
      <Q_QUAL_VALUE1>N</Q_QUAL_VALUE1>
    </ACTIVE>\n"""
        inactiveXml = etree.fromstring(inactiveStr)


        chemIdx = underlyingIdx.parent()
        dpIdx = underlyingIdx


        #add this to the appropriate element of the model. Need to
        #get the model index from the proxy model.
        
        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = chemIdx.internalPointer()
        dpti = dpIdx.internalPointer()

        #print 'chemti: ',chemti,chemti.parent()
        #print 'dpti: ',dpti

        sourceRoot = self.model.createIndex(0,0,self.model.rootItem).\
                child(0,0)
        proxyRoot = self.chemProxyModel.mapFromSource(sourceRoot)

        #for row in xrange(self.cv.chemProxyModel.rowCount(proxyRoot)):


        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()

        dpti.element.append(inactiveXml)
        self.model.treestep(inactiveXml,dpti)

        #print etree.tostring(self.model.rootItem.element,pretty_print=True)
        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()

        #indexes have all changed. Reset the root indices of the views.
        #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
        #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
        #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
        #self.chemView.setRootIndex(QtCore.QModelIndex())
        #self.dataPointView.setRootIndex(QtCore.QModelIndex().child(0,0))
        #self.qualifierView.setRootIndex(QtCore.QModelIndex().child(0,0).child(0,0))
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))


        #set the current index of the chemical view to the 
        #index of chemti. Need to make new index in self.model and
        #map to self.proxyModel.
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
        #similarly, set current index of data point view
        #self.dataPointSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
        #set current index of data point view to newly created
        #data point.
        self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)

        #print etree.tostring(xmlCopy,pretty_print=True)

        #get index of new datapoint

        #outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1]))

        outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti))


        #dummy - want to return the dataPointProxyModel index for 
        #the new data point.

        return outIdx

    def markAllAsInactive(self,idx):
        """Marks all datapoints of the same sort as the datapoint at 
        idx as inactive, for all visible chemicals."""
        

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)

        ti = underlyingIdx.internalPointer()
        xmlCopy = copy.deepcopy(ti.element)
        
        inactiveStr = """    <ACTIVE>
      <QUALIFIER_ID>AUTO</QUALIFIER_ID>
      <DATA_POINT_ID>AUTO</DATA_POINT_ID>
      <Q_L1_KINGDOM_CODE>Q_YN</Q_L1_KINGDOM_CODE>
      <CREATED_BY>AUTO</CREATED_BY>
      <CREATED_DATE>AUTO</CREATED_DATE>
      <Q_CONFIDENTIAL1_YN>N</Q_CONFIDENTIAL1_YN>
      <Q_QUAL_VALUE1>N</Q_QUAL_VALUE1>
    </ACTIVE>\n"""
        inactiveXml = etree.fromstring(inactiveStr)


        chemIdx = underlyingIdx.parent()
        dpIdx = underlyingIdx


        #add this to the appropriate element of the model. Need to
        #get the model index from the proxy model.
        
        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = chemIdx.internalPointer()
        dpti = dpIdx.internalPointer()

        #print 'chemti: ',chemti,chemti.parent()
        #print 'dpti: ',dpti

        sourceRoot = self.model.createIndex(0,0,self.model.rootItem).\
                child(0,0)
        proxyRoot = self.chemProxyModel.mapFromSource(sourceRoot)

        #for row in xrange(self.cv.chemProxyModel.rowCount(proxyRoot)):

        dataType = dpti.element.tag


        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()

        changedCount = 0

        #iterate over visible chemicals
        for row in xrange(self.chemProxyModel.rowCount(proxyRoot)):
            sourceVisTiIdx = proxyRoot.model().mapToSource(proxyRoot.child(row,0))
            visTi = sourceVisTiIdx.internalPointer()

            #Find all data points of the same type as the data point at idx
            for row2 in range(len(visTi.childItems)):
                if visTi.childItems[row2].element.tag == dataType:
                    #Is it active? If so, make it inactive.
                    if visTi.childItems[row2].element.xpath('./ACTIVE[last()]/Q_QUAL_VALUE1/text()=\'Y\''):
                        xmlCopy = copy.deepcopy(inactiveXml)
                        visTi.childItems[row2].element.append(xmlCopy)
                        self.model.treestep(xmlCopy,visTi.childItems[row2])

                        changedCount += 1
            

        #print etree.tostring(self.model.rootItem.element,pretty_print=True)
        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()

        #indexes have all changed. Reset the root indices of the views.
        #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
        #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
        #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
        #self.chemView.setRootIndex(QtCore.QModelIndex())
        #self.dataPointView.setRootIndex(QtCore.QModelIndex().child(0,0))
        #self.qualifierView.setRootIndex(QtCore.QModelIndex().child(0,0).child(0,0))
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))


        #set the current index of the chemical view to the 
        #index of chemti. Need to make new index in self.model and
        #map to self.proxyModel.
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
        #similarly, set current index of data point view
        #self.dataPointSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
        #set current index of data point view to newly created
        #data point.
        self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)

        #print etree.tostring(xmlCopy,pretty_print=True)

        #get index of new datapoint

        #outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1]))

        outIdx = self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti))


        #dummy - want to return the dataPointProxyModel index for 
        #the new data point.

        print changedCount+' data points de-activated. '

        return outIdx


    def editOrCopyQualifier(self,idx):
        """Edit qualifier or make a copy for editing. 
        
        If QUALIFIER_ID == "AUTO", really edit the qualifier, since it hasn't
        been uploaded to the remote database yet. Otherwise, make a copy of the 
        datapoint that points to the old one, so that there is always a record
        of the history of every piece of data."""
        #get DATA_POINT_ID. If it isn't 'AUTO', copy it. Otherwise,
        #edit it.

        #print 'qid: ', idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()
        qid = idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'QUALIFIER_ID')].toString()
        print 'qid: ',[qid]

        if qid == 'AUTO':
            self.editQualifier(idx)
        else:
            qName = str(idx.model().mapToSource(idx).data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'DP_L0_DOMAIN_CODE')].toPyObject())
            dpIdx = idx.model().mapToSource(idx).parent()
            if self.qualifierAllowed(qName,dpIdx):
                copyIdx = self.copyQualifier(idx)
                self.editQualifier(copyIdx)
   
    
    def editQualifier(self,idx):
        """Edit Qualifier. 
        
        Edit the qualifier, since it hasn't been uploaded to the 
        remote database yet"""

        print 'editQualifier ',idx

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)
        print underlyingIdx

        #make 'normal' dictionary (strings) from Qt types for input to
        #dataDialog.

        tempDict = {}
        for key in underlyingIdx.data(QtCore.Qt.UserRole).toMap().keys():
            tempDict[str(key)] = str(underlyingIdx.data(QtCore.Qt.UserRole).toMap()[key].toPyObject())


        #dg = dataDialog(self,underlyingIdx.data(QtCore.Qt.UserRole).toMap())
        dg = qualifierDialog(self,self.structureDict,tempDict)
        dgStatus = dg.exec_()
        if dgStatus:
            dgData = dg.getData()
            print '---'
            print dgData
            print '---'

            underlyingIdx.model().setData(underlyingIdx,dgData,QtCore.Qt.UserRole)

    def dataPointAllowed(self,itemName,chemIdx):
        """Check to see that adding a new datapoint of name itemName  is allowed."""
        print 'itemName',itemName
        rule = self.structureDict[itemName]['optional']
        #now, count how many of this sort of data point there are for this chemical.
        chemXml = chemIdx.internalPointer().element

        elementsWithSameName = chemXml.xpath('./'+itemName)
        no = len(elementsWithSameName)
        print rule,no

        if rule == 'one':
            if no == 1:
                return True
            else:
                return False
        elif rule == 'oneOrMore':
            if no >= 0:
                return True
            else:
                return False
        elif rule == 'zeroOrMore':
            return True
        elif rule == 'optional':
            if no == 0:
                return True
            else:
                return False
        else:
            #???
            return False


    def qualifierAllowed(self,qName,dpIdx):
        """Checks to see that adding a new qualifier of name qName  is allowed."""
        print 'qName',qName
        
        dpName = str(dpIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'DP_L0_DOMAIN_CODE')].toPyObject())
        
        rule = self.structureDict[dpName]['qualifiers'][qName]['optional']
        #now, count how many of this sort of data point there are for this chemical.
        dpXml = dpIdx.internalPointer().element

        elementsWithSameName = dpXml.xpath('./'+qName)
        no = len(elementsWithSameName)
        print rule,no

        if rule == 'one':
            if no == 1:
                return True
            else:
                return False
        elif rule == 'oneOrMore':
            if no >= 0:
                return True
            else:
                return False
        elif rule == 'zeroOrMore':
            return True
        elif rule == 'optional':
            if no == 0:
                return True
            else:
                return False
        else:
            #???
            return False

    def copyQualifier(self,idx):
        """Copy qualifier.
        
        Invoked when the user double clicks on a qualifier to edit it, and 
        the qualifier has been downloaded from the database. A copy of the 
        qualifier is created for the user to edit."""
        print 'copyQualifier'

        #get datapoint from underlying model
        underlyingIdx = idx.model().mapToSource(idx)
        print 'copyQualifier: ',underlyingIdx

        ti = underlyingIdx.internalPointer()
        xmlCopy = copy.deepcopy(ti.element)
        
        
        dpidXPath = etree.XPath('.//DATA_POINT_ID')
        dpids = dpidXPath(xmlCopy)

        for item in dpids:
            item.text = 'AUTO'


        qidXPath = etree.XPath('.//QUALIFIER_ID')
        qids = qidXPath(xmlCopy)

        for item in qids:
            item.text = 'AUTO'

        cbXPath = etree.XPath('.//CREATED_BY')
        cbs = cbXPath(xmlCopy)

        for item in cbs:
            item.text = 'AUTO'

        cdXPath = etree.XPath('.//CREATED_DATE')
        cds = cdXPath(xmlCopy)

        for item in cds:
            item.text = 'AUTO'
        print 'copyQualifier: ',etree.tostring(xmlCopy,pretty_print=True)


        chemIdx = underlyingIdx.parent().parent()
        dpIdx = underlyingIdx.parent()
        qIdx = underlyingIdx


        #add this to the appropriate element of the model. Need to
        #get the model index from the proxy model.
        
        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = chemIdx.internalPointer()
        dpti = dpIdx.internalPointer()
        qti = qIdx.internalPointer()


        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()
        dpti.element.append(xmlCopy)
        self.model.treestep(xmlCopy,dpti)

        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()

        #reset views
        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
        self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
        self.qualifierSelectionModel.setCurrentIndex(self.qualifierProxyModel.mapFromSource(self.model.createIndex(dpti.childItems[-1].row(),0,dpti.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)

        print etree.tostring(xmlCopy,pretty_print=True)

        #get index of new qualifier
        outIdx = self.qualifierProxyModel.mapFromSource(self.model.createIndex(dpti.childItems[-1].row(),0,dpti.childItems[-1]))


        #dummy - want to return the dataPointProxyModel index for 
        #the new data point.

        return outIdx


    
    def resizeEvent(self,event):
        """Emits resize signals."""
        self.model.emit(QtCore.SIGNAL('dataChanged'))
        self.model.emit(QtCore.SIGNAL('layoutChanged'))

    def makeContextMenus(self):
        """Makes datapoint and qualifier context menus."""
        self.chemContextMenu = QtGui.QMenu(self)
        self.chemContextMenu.addAction('Delete selected chemical',self.onActionDeleteChemical)
        self.dpContextMenu = QtGui.QMenu(self)
        self.dpContextMenu.addAction('Delete selected data point',self.onActionDeleteDp)

        self.chemView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.chemView.customContextMenuRequested.connect(self.onChemContextMenu)

        self.dataPointView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.dataPointView.customContextMenuRequested.connect(self.onDpContextMenu)

        self.qContextMenu = QtGui.QMenu(self)
        self.qContextMenu.addAction('Delete selected qualifier',self.onActionDeleteQ)
        self.qualifierView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.qualifierView.customContextMenuRequested.connect(self.onQContextMenu)
        print 'finished makeContextMenus'
    
    def onActionDeleteChemical(self):
        """Deletes chemical"""
        print 'delete chemical!'
        
        chemIdx = self.chemSelectionModel.currentIndex()
        chemIdx = self.chemProxyModel.mapToSource(chemIdx)

        dpIdx = self.dataPointSelectionModel.currentIndex()
        dpIdx = self.dataPointProxyModel.mapToSource(dpIdx)

        #Get selected chemical's xml element
        chemelement = chemIdx.internalPointer().element

        print chemelement

        #double check to make sure that DATA_POINT_ID == 'AUTO'
        CHEMICAL_ID = chemIdx.internalPointer().data(QtCore.Qt.UserRole)['CHEMICAL_ID']
        
        if CHEMICAL_ID == 'AUTO':
            print 'CHEMICAL_ID == AUTO'


            #save pointers to currently selected treeItems to reset later.
            chemti = chemIdx.internalPointer()
            chemtiChildNumber = chemIdx.row()
            dpti = dpIdx.internalPointer()
            dptiChildNumber = dpIdx.row()

            #DP_L0_DOMAIN_CODE=dpIdx.internalPointer().data(QtCore.Qt.UserRole)[QtCore.QString(u'DP_L0_DOMAIN_CODE')].toString()
            

            #work out where to set the selection to after
            #resetting the model.
            lastIdx = chemIdx.internalPointer().parent().childCount()
            lastIdx = lastIdx - 1

            if chemtiChildNumber != lastIdx:
                #The selected chemical is not the last in the list.
                #Use the next ti in the list
                finalSelectionTi = chemIdx.internalPointer().parent().childItems[chemtiChildNumber+1]
                finalSelectionIdx = chemtiChildNumber
            elif chemtiChildNumber == lastIdx and lastIdx > 0:
                #The selected chemical is the last in the list and 
                #there are at least two elements in the list.
                #Use the previous ti in the list.
                finalSelectionTi = chemIdx.internalPointer().parent().childItems[chemtiChildNumber-1]
                finalSelectionIdx = chemtiChildNumber-1
            else:
                pass

            self.model.beginResetModel()
            self.chemProxyModel.beginResetModel()
            self.dataPointProxyModel.beginResetModel()
            self.qualifierProxyModel.beginResetModel()

            chemelement.getparent().remove(chemelement)

            #Get index in childItems list of the parent of the treeitem to be
            #removed.
            childIdx = chemIdx.internalPointer().parent().childItems.index(chemIdx.internalPointer())
            del  chemIdx.internalPointer().parent().childItems[childIdx]

            self.model.endResetModel()
            self.chemProxyModel.endResetModel()
            self.dataPointProxyModel.endResetModel()
            self.qualifierProxyModel.endResetModel()
 
            #reset views.
            self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
            if lastIdx != 0:
                self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(finalSelectionTi.row(),0,finalSelectionTi)),QtGui.QItemSelectionModel.SelectCurrent)
            else:
                self.chemSelectionModel.clear()
                self.dataPointSelectionModel.clear()
                self.qualifierSelectionModel.clear()
                #Hack to make sure nothing is visible in any of the views.
                #Set the root indices of all of the views to one above the chemical
                #level.
                self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
                self.dataPointView.setRootIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
                self.qualifierView.setRootIndex(self.qualifierProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
            #refresh figure.
            #if DP_L0_DOMAIN_CODE == 'SMILES':
            #    self.model.makeOneSvg(chemti)
            #    #self.onSplitterMoved(1,QtCore.QModelIndex())
            #    self.chemViewDe.sizeHintChanged.emit(dpIdx.parent())



    def onActionDeleteDp(self):
        """Deletes data point"""
        print 'delete dp!'
        
        chemIdx = self.chemSelectionModel.currentIndex()
        chemIdx = self.chemProxyModel.mapToSource(chemIdx)

        dpIdx = self.dataPointSelectionModel.currentIndex()
        dpIdx = self.dataPointProxyModel.mapToSource(dpIdx)

        #Get selected data point xml element
        dpelement = dpIdx.internalPointer().element

        print dpelement

        #double check to make sure that DATA_POINT_ID == 'AUTO'
        DATA_POINT_ID = dpIdx.internalPointer().data(QtCore.Qt.UserRole)['DATA_POINT_ID']
        
        if DATA_POINT_ID == 'AUTO':


            #save pointers to currently selected treeItems to reset later.
            chemti = chemIdx.internalPointer()
            chemtiChildNumber = chemIdx.row()
            dpti = dpIdx.internalPointer()
            dptiChildNumber = dpIdx.row()

            DP_L0_DOMAIN_CODE=dpIdx.internalPointer().data(QtCore.Qt.UserRole)['DP_L0_DOMAIN_CODE']
            

            self.model.beginResetModel()
            self.chemProxyModel.beginResetModel()
            self.dataPointProxyModel.beginResetModel()
            self.qualifierProxyModel.beginResetModel()
            
            dpelement.getparent().remove(dpelement)

            #Get index in childItems list of the parent of the treeitem to be
            #removed.
            childIdx = dpIdx.internalPointer().parent().childItems.index(dpIdx.internalPointer())
            del  dpIdx.internalPointer().parent().childItems[childIdx]

            self.model.endResetModel()
            self.chemProxyModel.endResetModel()
            self.dataPointProxyModel.endResetModel()
            self.qualifierProxyModel.endResetModel()
 
            #reset views.
            self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))

            self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
            #refresh figure.
            if DP_L0_DOMAIN_CODE == 'SMILES':
                self.model.makeOneSvg(chemti)
                #self.onSplitterMoved(1,QtCore.QModelIndex())
                self.chemViewDe.sizeHintChanged.emit(dpIdx.parent())



    def onActionDeleteQ(self):
        """Deletes qualifier"""
        print 'delete qualifier!'
        
        chemIdx = self.chemSelectionModel.currentIndex()
        chemIdx = self.chemProxyModel.mapToSource(chemIdx)

        dpIdx = self.dataPointSelectionModel.currentIndex()
        dpIdx = self.dataPointProxyModel.mapToSource(dpIdx)
        
        qIdx = self.qualifierSelectionModel.currentIndex()
        qIdx = self.qualifierProxyModel.mapToSource(qIdx)

        #Get selected qualifier xml elementa
        print qIdx.internalPointer().data(QtCore.Qt.UserRole)
        qelement = qIdx.internalPointer().element

        print qelement
        dpTag = str(qelement.getparent().tag)
        qTag = str(qelement.tag)

        #double check to make sure that DATA_POINT_ID == 'AUTO'
        QUALIFIER_ID = qIdx.internalPointer().data(QtCore.Qt.UserRole)['QUALIFIER_ID']
        
        if QUALIFIER_ID == 'AUTO':


            #save pointers to currently selected treeItems to reset later.
            chemti = chemIdx.internalPointer()
            chemtiChildNumber = chemIdx.row()
            dpti = dpIdx.internalPointer()
            dptiChildNumber = dpIdx.row()
            qti = qIdx.internalPointer()
            qtiChildNumber = qIdx.row()

            self.model.beginResetModel()
            self.chemProxyModel.beginResetModel()
            self.dataPointProxyModel.beginResetModel()
            self.qualifierProxyModel.beginResetModel()

            qelement.getparent().remove(qelement)


            #Get index in childItems list of the parent of the treeitem to be
            #removed.
            childIdx = qIdx.internalPointer().parent().childItems.index(qIdx.internalPointer())
            del  qIdx.internalPointer().parent().childItems[childIdx]

            self.model.endResetModel()
            self.chemProxyModel.endResetModel()
            self.dataPointProxyModel.endResetModel()
            self.qualifierProxyModel.endResetModel()
 
            #reset views.
            self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))

            self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)

            self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
            if dpTag == 'SMILES' and \
                    qTag == 'ACTIVE':
                self.model.makeOneSvg(chemIdx.internalPointer())
                self.onSplitterMoved(1,QtCore.QModelIndex())
                self.chemViewDe.sizeHintChanged.emit(chemIdx)


    def onDpContextMenu(self,point):
        """Slot which makes a data point context menu."""

        #get the data point under the mouse pointer.
        a = self.dataPointProxyModel.mapToSource(self.dataPointView.indexAt(point))
        
        if a.internalPointer() != None:
            #check to see that DATA_POINT_ID == 'AUTO'
            DATA_POINT_ID = a.internalPointer().data(QtCore.Qt.UserRole)['DATA_POINT_ID']
            print DATA_POINT_ID

            if DATA_POINT_ID == 'AUTO':
                self.dpContextMenu.exec_(self.dataPointView.mapToGlobal(point))

    def onChemContextMenu(self,point):
        """Slot which makes a chemical context menu."""

        #get the data point under the mouse pointer.
        a = self.chemProxyModel.mapToSource(self.chemView.indexAt(point))
        
        if a.internalPointer() != None:
            #check to see that DATA_POINT_ID == 'AUTO'
            CHEMICAL_ID = a.internalPointer().data(QtCore.Qt.UserRole)['CHEMICAL_ID']
            print CHEMICAL_ID

            if CHEMICAL_ID == 'AUTO':
                self.chemContextMenu.exec_(self.chemView.mapToGlobal(point))


    def onQContextMenu(self,point):
        """Slot which makes qualifier context menu."""

        #get the data point under the mouse pointer.
        a = self.qualifierProxyModel.mapToSource(self.qualifierView.indexAt(point))
        
        if a.internalPointer() != None:
            #check to see that QUALIFIER_ID == 'AUTO'
            QUALIFIER_ID = a.internalPointer().data(QtCore.Qt.UserRole)['QUALIFIER_ID']
            print QUALIFIER_ID

            if QUALIFIER_ID == 'AUTO':
                self.qContextMenu.exec_(self.qualifierView.mapToGlobal(point))

    
    def newChemical(self):
        """Adds a new chemical."""
        print 'New Chemical'
        newxml = etree.Element('CHEMICAL')
        chemicalId = etree.SubElement(newxml,'CHEMICAL_ID')
        chemicalId.text = 'AUTO'
        createdBy = etree.SubElement(newxml,'CREATED_BY')
        createdBy.text = 'AUTO'
        createdDate = etree.SubElement(newxml,'CREATED_DATE')
        createdDate.text = 'AUTO'
        print 'new xml created'

        chemIdx = self.chemSelectionModel.currentIndex()
        dpIdx = self.dataPointSelectionModel.currentIndex()
        qIdx = self.qualifierSelectionModel.currentIndex()
        print 'saved indices'

        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = self.chemProxyModel.mapToSource(chemIdx).internalPointer()
        dpti = self.dataPointProxyModel.mapToSource(dpIdx).internalPointer()
        qti = self.qualifierProxyModel.mapToSource(qIdx).internalPointer()

        print 'saved treeitems'

        self.model.beginResetModel()
        self.chemProxyModel.beginResetModel()
        self.dataPointProxyModel.beginResetModel()
        self.qualifierProxyModel.beginResetModel()


        #self.chemProxyModel.mapToSource(chemIdx).parent().internalPointer().element.append(newxml)
        print 'newChemical: self.model.rootItem.element',self.model.rootItem.element
        print 'newChemical: self.model.rootItem.childCount()',self.model.rootItem.childCount()

        #hack to add a root item if it doesn't exist. treestep doesn't work
        #if the <ROOT> element doesn't contain any children - it thinks the <ROOT>
        #element is a leaf rather than a real element.
        if self.model.rootItem.childCount() == 0:
            rootEl = self.model.rootItem.element.xpath('ROOT')[0]
            print rootEl
            newitem = TreeItem(rootEl,self.model.rootItem)
            self.model.rootItem.appendChild(newitem)

        self.model.rootItem.child(0).element.append(newxml)
        self.model.treestep(newxml,self.model.rootItem.child(0))
        
        #print etree.tostring(self.model.rootItem.element,pretty_print=True)
        self.model.endResetModel()
        self.chemProxyModel.endResetModel()
        self.dataPointProxyModel.endResetModel()
        self.qualifierProxyModel.endResetModel()


        self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
        print 'new xml added.'
        print 'foo'
        newTi = self.model.rootItem.child(0).childItems[-1]
        self.model.makeOneSvg(newTi)
        print newTi.data(0)
        #print etree.tostring(newTi.parent().data(0)[QtCore.QString(u'xmlelement')].toPyObject(),pretty_print=True)
        self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(newTi.row(),0,newTi)),QtGui.QItemSelectionModel.SelectCurrent)
        #indexes have all changed. Reset the root indices of the views.
        #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
        #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
        #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
        #
        #print '>>3'
        #self.chemSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(self.model.rootItem.childItems[-1].row(),0,self.model.rootItem.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)
        #print '>>4'

        newIdx = self.model.createIndex(newTi.row(),0,newTi)
        self.model.dataChanged.emit(newIdx,newIdx)



    def newDataPoint(self):
        """Add a new datapoint.
        
        Makes a dataDialog prompting the user to enter a new piece of data. The 
        data collected by the dialog is then added to the underlying xml, and
        the model is updated using TreeModel.treestep.
        
        Note that this duplicates functionality found in TreeModel.setData."""
        print 'New Data Point'
        chemIdx = self.chemSelectionModel.currentIndex()
        dpIdx = self.dataPointSelectionModel.currentIndex()
        qIdx = self.qualifierSelectionModel.currentIndex()

        #firstly, save pointers to TreeItems so that the views
        #can be reset
        chemti = self.chemProxyModel.mapToSource(chemIdx).internalPointer()
        print 'chemti'
        dpti = self.dataPointProxyModel.mapToSource(dpIdx).internalPointer()
        print dpIdx
        print 'dpti'
        qti = self.dataPointProxyModel.mapToSource(qIdx).internalPointer()
        print 'qti'

        dg = dataDialog(self,self.structureDict)
        dgStatus = dg.exec_()
        if dgStatus:
            dgData = dg.getData()
            print dgData
            #make piece of xml to tack onto what already exists.
            newxml = etree.Element(dgData['DP_L0_DOMAIN_CODE'])
            #print etree.tostring(newxml,pretty_print=True)
            
            DATA_POINT_ID = etree.SubElement(newxml,'DATA_POINT_ID')
            DATA_POINT_ID.text = 'AUTO'

            #get index of chemical currently selected in self.chemView
            #This index refers to proxyModel. Get data as a dict, and 
            #add the text to the CHEMICAL_ID element. Note that chemIdx
            #is for the proxy model.
            #chemIdx = self.chemView.currentIndex()
            #dpIdx = self.dataPointView.currentIndex()

            #print 'dpIdx: ',self.proxyModel.mapToSource(dpIdx).internalPointer()
            #print 'chemIdx: ',chemIdx
            #print 'qIdx: ',qIdx
            CHEMICAL_ID = etree.SubElement(newxml,'CHEMICAL_ID')
            print chemIdx.data().toMap()
            CHEMICAL_ID.text = str(chemIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('CHEMICAL_ID')].toPyObject())

            print '----'

            if 'DP_L1_KINGDOM_CODE' in dgData.keys():
                DP_L1_KINGDOM_CODE = etree.SubElement(newxml,'DP_L1_KINGDOM_CODE')
                DP_L1_KINGDOM_CODE.text = dgData['DP_L1_KINGDOM_CODE']
            else:
                print 'DP_L1_KINGDOM_CODE not in dgData'

            
            CREATED_BY = etree.SubElement(newxml,'CREATED_BY')
            CREATED_BY.text = 'AUTO'
            
            CREATED_DATE = etree.SubElement(newxml,'CREATED_DATE')
            CREATED_DATE.text = 'AUTO'
            
            if 'DP_QUANT_VALUE1' in dgData.keys():
                DP_QUANT_VALUE1 = etree.SubElement(newxml,'DP_QUANT_VALUE1')
                DP_QUANT_VALUE1.text = dgData['DP_QUANT_VALUE1']
            else:
                print 'DP_QUANT_VALUE1 not in dgData'
            
            if 'DP_QUANT_VALUE2' in dgData.keys():
                DP_QUANT_VALUE2 = etree.SubElement(newxml,'DP_QUANT_VALUE2')
                DP_QUANT_VALUE2.text = dgData['DP_QUANT_VALUE2']
            else:
                print 'DP_QUANT_VALUE2 not in dgData'


            if 'DP_INTEGER1' in dgData.keys():
                DP_INTEGER1 = etree.SubElement(newxml,'DP_INTEGER1')
                DP_INTEGER1.text = dgData['DP_INTEGER1']
            else:
                print 'DP_INTEGER1 not in dgData'
            
            if 'DP_INTEGER2' in dgData.keys():
                DP_INTEGER2 = etree.SubElement(newxml,'DP_INTEGER2')
                DP_INTEGER2.text = dgData['DP_INTEGER2']
            else:
                print 'DP_INTEGER1 not in dgData'
            
            if 'DP_CLOB1' in dgData.keys():
                DP_CLOB1 = etree.SubElement(newxml,'DP_CLOB1')
                DP_CLOB1.text = dgData['DP_CLOB1']
            else:
                print 'DP_CLOB1 not in dgData'
            
            if 'DP_CONFIDENTIAL1_YN' in dgData.keys():
                DP_INTEGER2 = etree.SubElement(newxml,'DP_CONFIDENTIAL1_YN')
                DP_INTEGER2.text = dgData['DP_CONFIDENTIAL1_YN']
            else:
                print 'DP_CONFIDENTIAL1_YN not in dgData'
            
            if 'DP_CONFIDENTIAL2_YN' in dgData.keys():
                DP_INTEGER2 = etree.SubElement(newxml,'DP_CONFIDENTIAL2_YN')
                DP_INTEGER2.text = dgData['DP_CONFIDENTIAL2_YN']
            else:
                print 'DP_CONFIDENTIAL2_YN not in dgData'

            if 'DP_QUAL_VALUE1' in dgData.keys():
                DP_QUAL_VALUE1 = etree.SubElement(newxml,'DP_QUAL_VALUE1')
                DP_QUAL_VALUE1.text = dgData['DP_QUAL_VALUE1']
            else:
                print 'DP_QUAL_VALUE1 not in dgData'
            
            if 'DP_QUAL_VALUE2' in dgData.keys():
                DP_QUAL_VALUE2 = etree.SubElement(newxml,'DP_QUAL_VALUE2')
                DP_QUAL_VALUE2.text = dgData['DP_QUAL_VALUE2']
            else:
                print 'DP_QUAL_VALUE2 not in dgData'

            if 'DP_UNIT_CODE1' in dgData.keys():
                DP_UNIT_CODE1 = etree.SubElement(newxml,'DP_UNIT_CODE1')
                DP_UNIT_CODE1.text = dgData['DP_UNIT_CODE1']
            else:
                print 'DP_UNIT_CODE2 not in dgData'

            if 'DP_UNIT_CODE2' in dgData.keys():
                DP_UNIT_CODE2 = etree.SubElement(newxml,'DP_UNIT_CODE2')
                DP_UNIT_CODE2.text = dgData['DP_UNIT_CODE2']
            else:
                print 'DP_UNIT_CODE2 not in dgData'
            

            #Add ACTIVE qualifier (datapoints are active by default when added).

            ACTIVE = etree.SubElement(newxml,'ACTIVE')
            qidAuto = etree.SubElement(ACTIVE,'QUALIFIER_ID')
            qidAuto.text = 'AUTO'
            qidAuto = etree.SubElement(ACTIVE,'DATA_POINT_ID')
            qidAuto.text = 'AUTO'
            ql1kc = etree.SubElement(ACTIVE,'Q_L1_KINGDOM_CODE')
            ql1kc.text ='Q_YN'
            qcb = etree.SubElement(ACTIVE,'CREATED_BY')
            qcb.text = 'AUTO'
            qcd = etree.SubElement(ACTIVE,'CREATED_DATE')
            qcd.text = 'AUTO'
            qconf = etree.SubElement(ACTIVE,'Q_CONFIDENTIAL1_YN')
            qconf.text = 'N'
            qqv1 = etree.SubElement(ACTIVE,'Q_QUAL_VALUE1')
            qqv1.text = 'Y'

            if dgData['DP_L0_DOMAIN_CODE'] == 'Endpoint':
                #set default pivotal endpoint tag to 'N'
                PIVOTAL = etree.SubElement(newxml,'PIVOTAL')
                qidAuto = etree.SubElement(PIVOTAL,'QUALIFIER_ID')
                qidAuto.text = 'AUTO'
                qidAuto = etree.SubElement(PIVOTAL,'DATA_POINT_ID')
                qidAuto.text = 'AUTO'
                ql1kc = etree.SubElement(PIVOTAL,'Q_L1_KINGDOM_CODE')
                ql1kc.text ='Q_YN'
                qcb = etree.SubElement(PIVOTAL,'CREATED_BY')
                qcb.text = 'AUTO'
                qcd = etree.SubElement(PIVOTAL,'CREATED_DATE')
                qcd.text = 'AUTO'
                qconf = etree.SubElement(PIVOTAL,'Q_CONFIDENTIAL1_YN')
                qconf.text = 'N'
                qqv1 = etree.SubElement(PIVOTAL,'Q_QUAL_VALUE1')
                qqv1.text = 'N'

            print '----'

            #add this to the appropriate element of the model. Need to
            #get the model index from the proxy model.
            

            #print 'chemti: ',chemti,chemti.parent()
            #print 'dpti: ',dpti

            self.model.beginResetModel()
            self.chemProxyModel.beginResetModel()
            self.dataPointProxyModel.beginResetModel()
            self.qualifierProxyModel.beginResetModel()

            self.chemProxyModel.mapToSource(chemIdx).internalPointer().element.append(newxml)
            self.model.treestep(newxml,self.chemProxyModel.mapToSource(chemIdx).internalPointer())
            #print etree.tostring(self.model.rootItem.element,pretty_print=True)
            self.model.endResetModel()
            self.chemProxyModel.endResetModel()
            self.dataPointProxyModel.endResetModel()
            self.qualifierProxyModel.endResetModel()

            #indexes have all changed. Reset the root indices of the views.
            #self.chemView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()))
            #self.dataPointView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0))
            #self.qualifierView.setRootIndex(self.proxyModel.index(0,0,QtCore.QModelIndex()).child(0,0).child(0,0))
            #self.chemView.setRootIndex(QtCore.QModelIndex())
            #self.dataPointView.setRootIndex(QtCore.QModelIndex().child(0,0))
            #self.qualifierView.setRootIndex(QtCore.QModelIndex().child(0,0).child(0,0))
            self.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))

            #set the current index of the chemical view to the 
            #index of chemti. Need to make new index in self.model and
            #map to self.proxyModel.
            self.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
            #similarly, set current index of data point view
            #self.dataPointSelectionModel.setCurrentIndex(self.proxyModel.mapFromSource(self.model.createIndex(dpti.row(),0,dpti)),QtGui.QItemSelectionModel.SelectCurrent)
            #set current index of data point view to newly created
            #data point.
            self.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)

            print etree.tostring(newxml,pretty_print=True)

            if dgData['DP_L0_DOMAIN_CODE'] == 'SMILES':
                self.model.makeOneSvg(chemti)
                self.onSplitterMoved(1,QtCore.QModelIndex())
                self.chemViewDe.sizeHintChanged.emit(chemIdx)

            newDataPointIndex = self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1])

            self.model.dataChanged.emit(newDataPointIndex,newDataPointIndex)



    
    
    def newQualifier(self):
        """Adds a new qualifier.
        
        Note that this duplicates functionality in TreeModel.setData."""
        def insertXml(dp,name,inDict,insertPosition):
            if name in inDict.keys():
                a = etree.Element(name)
                a.text = inDict[name]
                dp.insert(insertPosition,a)
                insertPosition += 1
            return insertPosition


        print 'New Qualifier'

        dpIdx = self.dataPointProxyModel.mapToSource(\
                self.dataPointSelectionModel.currentIndex())
        dpti = dpIdx.internalPointer()

        qIdx = self.qualifierProxyModel.mapToSource(\
                self.qualifierSelectionModel.currentIndex())
        qti = qIdx.internalPointer()

        q = qualifierDialog(self,self.structureDict)
        qStatus = q.exec_()
        if qStatus:
            qData = q.getData()
            print qData
            
            #Note that 'DP_L0_DOMAIN_CODE' is used for qualifiers 
            #as well, for historical reasons.
            newxml = etree.Element(qData['DP_L0_DOMAIN_CODE'])
            #print etree.tostring(newxml,pretty_print=True)
            
            QUALIFIER_ID = etree.SubElement(newxml,'QUALIFIER_ID')
            QUALIFIER_ID.text = 'AUTO'

            DATA_POINT_ID = etree.SubElement(newxml,'DATA_POINT_ID')
            DATA_POINT_ID.text = str(dpIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString('DATA_POINT_ID')].toPyObject())


            insertPosition = insertXml(newxml,'Q_L1_KINGDOM_CODE',qData,2)

            CREATED_BY = etree.Element('CREATED_BY')
            CREATED_BY.text = 'AUTO'
            newxml.insert(insertPosition,CREATED_BY)
            insertPosition += 1

            CREATED_DATE = etree.Element('CREATED_DATE')
            CREATED_DATE.text = 'AUTO'
            newxml.insert(insertPosition,CREATED_DATE)
            insertPosition += 1

            eList = ['Q_QUANT_VALUE1','Q_QUANT_VALUE2',\
                'Q_CLOB1', 'Q_INTEGER1',\
                'Q_CONFIDENTIAL1_YN', 'Q_CONFIDENTIAL2_YN',\
                'Q_QUAL_VALUE1','Q_QUAL_VALUE2',\
                'Q_UNIT_CODE1','Q_UNIT_CODE2']

            for name in eList:
                insertPosition = insertXml(newxml,name,qData,insertPosition)

            if dpti.childCount() == 0:
                newRow = 0
            elif dpti.childCount() > 0:
                newRow = qIdx.parent().internalPointer().childCount()
            else:
                print 'newRow: ???'

            dpti.element.append(newxml)

            self.model.beginInsertRows(dpIdx,newRow,newRow)
            self.model.treestep(newxml,dpti)
            self.model.endInsertRows()
            newIdx = self.model.index(newRow,0,dpIdx)
            self.model.dataChanged.emit(newIdx,newIdx)

            #print etree.tostring(dpIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'xmlelement')].toPyObject(),pretty_print=True)




            self.qualifierSelectionModel.clearSelection()
            self.qualifierSelectionModel.setCurrentIndex(self.qualifierProxyModel.mapFromSource(self.model.createIndex(dpti.childItems[-1].row(),0,dpti.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)

            idxOfNewQualifier = self.model.createIndex(dpti.childItems[-1].row(),0,dpti.childItems[-1])
            self.model.dataChanged.emit(idxOfNewQualifier,idxOfNewQualifier)
            
            if str(newxml.tag) == 'ACTIVE' and \
                    str(newxml.getparent().tag) == 'SMILES':
                self.model.makeOneSvg(dpti.parent())
                #self.onSplitterMoved(1,QtCore.QModelIndex())
                self.chemViewDe.sizeHintChanged.emit(dpIdx.parent())
            #qIdx.internalPointer().data(QtCore.Qt.UserRole)[QtCore.QString(u'xmlelement')].toPyObject()

#

    def dictFromSchema(self):
        """Gets a dictionary specifying xml structure from relaxNG schema."""
        
        print 'dictFromSchema'
        
        #get datapoint schema entries.
        dps = self.schemaXML.xpath('//rng:start/rng:element[@name=\'ROOT\']'\
                '/rng:zeroOrMore/rng:element[@name=\'CHEMICAL\']/rng:interleave'\
                '/rng:zeroOrMore/rng:element[@name]',\
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})

        dpDict = {}
        #get names
        for dp in dps:
            #qualifiers: subelements of dp that have a <ref name="QHEADER">
            #subelement.
            dpDict[dp.get('name')] = {}

            #get long name from comment field
            longName = dp.xpath('./comment()')
            #print longName[0].text
            dpDict[dp.get('name')]['longName'] = longName[0].text

            #work out whether it's optional or not
            dpZeroOrMore = dp.xpath('boolean(./parent::rng:zeroOrMore)',
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
            dpOptional = dp.xpath('boolean(./parent::rng:optional)',
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
            dpOneOrMore =dp.xpath('boolean(./parent::rng:oneOrMore)',
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})

            if dpZeroOrMore == True:
                dpDict[dp.get('name')]['optional'] = 'zeroOrMore'
            elif dpOptional == True:
                 dpDict[dp.get('name')]['optional'] = 'optional'
            elif dpOneOrMore == True:
                 dpDict[dp.get('name')]['optional'] = 'oneOrMore'
            else: # it's not optional
                dpDict[dp.get('name')]['optional'] = 'one'

            #get dpType
            regExpNs = 'http://exslt.org/regular-expressions'
            dpTypes = dp.xpath('.//rng:ref[re:test(./@name,\'^DP_.*$\')]',
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                're' :regExpNs})
            
            if len(dpTypes) > 1:
                floatSet = set([])
                for item in dpTypes:
                    floatSet.add(item.get('name'))
                #print floatSet
                if floatSet == set(['DP_FLOAT_GT', 'DP_FLOAT_SINGLE', 'DP_FLOAT_RANGE', 'DP_FLOAT_LT', 'DP_FLOAT_EXCEPTION']) or floatSet ==set(['DP_FLOAT_GT', 'DP_FLOAT_SINGLE', 'DP_FLOAT_RANGE', 'DP_FLOAT_LT']):
                    dpType = 'DP_FLOAT'
                else:
                    print 'len(dpTypes) > 1 and not a float' 
                    print floatSet
                    exit()
                
                if 'DP_FLOAT_EXCEPTION' in floatSet:
                    #get the values of DP_VARCHAR2 that are allowed
                    dpExcept = dp.xpath('.//rng:element[@name=\'DP_QUAL_VALUE1\']/rng:value',
                            namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'\
                            , 're' :regExpNs})

                    exceptList = []
                    for item in dpExcept:
                        exceptList.append(item.text)

                    dpDict[dp.get('name')]['exceptions'] = exceptList

           
            elif len(dpTypes) == 1:
                dpType = dpTypes[0].get('name')
            
            else:
                print 'len(dpTypes) != 1'
                print 'dpTypes: ', dpTypes
                sys.exit()

            dpDict[dp.get('name')]['dpType'] = dpType
            
            #If it's a float, get allowed units.
            if dpType == 'DP_FLOAT':
                #get units
                unit1 = dp.xpath('.//rng:element[@name=\'DP_UNIT_CODE1\']/rng:value|.//rng:element[@name=\'DP_UNIT_CODE1\']/rng:choice/rng:value',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
                if len(unit1) > 0:
                    unit1 = set([x.text for x in unit1])
                else:
                    unit1 = set([])
                #print 'unit1: ',unit1
                
                unit2 = dp.xpath('.//rng:element[@name=\'DP_UNIT_CODE2\']/rng:value|.//rng:element[@name=\'DP_UNIT_CODE2\']/rng:choice/rng:value',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
                if len(unit2) > 0:
                    unit2 = set([x.text for x in unit2])
                else:
                    unit2 = set([])
                #print 'unit2:', unit2
                
                if unit2 != unit1:
                    print 'unit 1 != unit2'
                    sys.exit()
                dpDict[dp.get('name')]['dpUnits'] = unit1
            elif dpType == 'DP_VARCHAR2CHOICE':
                #get choices
                choices=dp.xpath('.//rng:element[@name=\'DP_QUAL_VALUE1\']/rng:choice/rng:value/text()',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
                dpDict[dp.get('name')]['choices'] = choices
                choices.sort(key=lambda x: x.lower())
 
            
            #now get list of qualifier names to add to dpDict and a list of 
            #qualifiers to add to the qualifier dict
            #qualifiers = dp.xpath('.//rng:element[.//rng:ref[@name=\'QHEADER\']]',
            regExpNs = 'http://exslt.org/regular-expressions'
            qualifiers = dp.xpath('.//rng:element[descendant::rng:ref[re:test(./@name,\'^Q_.*$\')]]',\
                namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                're' :regExpNs})


            qualDict = {}
            dpDict[dp.get('name')]['qualifiers']={}
            for qualifier in qualifiers:
                #print qualifier
                dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')] = {}
                
                #work out whether it's optional or not
                qZeroOrMore = qualifier.xpath('boolean(./parent::rng:zeroOrMore)',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
                qOptional = qualifier.xpath('boolean(./parent::rng:optional)',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})
                qOneOrMore =qualifier.xpath('boolean(./parent::rng:oneOrMore)',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0'})

                if qZeroOrMore == True:
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['optional'] = 'zeroOrMore'
                elif qOptional == True:
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['optional'] =  'optional'
                elif qOneOrMore == True:
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['optional'] = 'oneOrMore'
                else: # it's not optional
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['optional'] = 'one'




                #get qType
                qTypes = qualifier.xpath('.//rng:ref[re:test(./@name,\'^Q_.*$\')]',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                    're' :regExpNs})


                if len(qTypes) > 1:
                    floatSet = set([])
                    for item in qTypes:
                        floatSet.add(item.get('name'))
                    #print floatSet
                    if floatSet == set(['Q_FLOAT_GT', 'Q_FLOAT_SINGLE', 'Q_FLOAT_RANGE', 'Q_FLOAT_LT']):
                        qType = 'Q_FLOAT'

                        #get units for float.
                        print etree.tostring(qualifier)
                        unit1 = qualifier.xpath('.//rng:element[@name="Q_UNIT_CODE1"]//rng:value',
                            namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                            're' :regExpNs})
                        unit1 = set([x.text for x in unit1])
                        unit2 = qualifier.xpath('.//rng:element[@name="Q_UNIT_CODE2"]//rng:value',
                            namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                            're' :regExpNs})
                        unit2 = set([x.text for x in unit2])
                        if unit1 != unit2:
                            print 'unit1 != unit2!'
                            print 'unit1:',unit1
                            print 'unit2:',unit2
                            exit()

                        dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                                ['qUnits'] = unit1




                    else:
                        'len(qTypes) > 1 and not a float'
                        exit()


                elif len(qTypes) == 1:
                    qType = qTypes[0].get('name')

                else:
                    print 'len(qTypes) != 1'
                    print 'qTypes: ', qTypes
                    sys.exit()

                dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['qType'] = qType




                #get long name from comment
                ln = qualifier.xpath('./comment()',
                    namespaces={'rng': 'http://relaxng.org/ns/structure/1.0',\
                    're' :regExpNs})
                if len(ln) >= 1:
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['longName'] = ln[0].text
                else:
                    dpDict[dp.get('name')]['qualifiers'][qualifier.get('name')]\
                        ['longName'] = None

        #pprint(dpDict)
        self.structureDict = dpDict
        #return dpDict




class dataDialog(QtGui.QDialog):
    """Dialog box to allow addition of a new data point.
    
    This is also used to make a copy of an old data point
    for further editing."""
    def __init__(self,parent,dpDict,initDict=None):
        super(dataDialog,self).__init__(parent)
        if initDict !=None:
            self.initDict = initDict
        else:
            self.initDict = None

        self.dpdict = dpDict

        print 'dataDialog.initDict: ', self.initDict

        #self.schemaXML = etree.parse('newschema.rng').getroot()
        #self.schemaRNG = etree.RelaxNG(self.schemaXML)
        
        self.parentChemIdx = parent.chemProxyModel.mapToSource(parent.chemSelectionModel.currentIndex())

        self.initUI()

    def initUI(self):
        """Initialises user interface of dialog."""

        #make combo box.
        self.dataPointCombo = QtGui.QComboBox()
        #items = ['DP_YN','DP_YNPStar','DP_CLOB','DP_FLOAT']
        #self.dpdict = self.dictFromSchema()
        #pprint(self.dpdict)
        #print 'Done dpDict'
       
        #list all possible entries
        items = self.dpdict.keys()
        items.sort(key=str.lower)
        
        #count all of the datapoints the chemical currently has.
        itemCount = {}
        for item in self.parentChemIdx.internalPointer().childItems:
            itemName = item.datadict['DP_L0_DOMAIN_CODE']
            if itemName not in itemCount.keys():
                itemCount[itemName] = 1
            else:
                itemCount[itemName] += 1

        #For everything else in items, set itemCount to zero.
        for item in items:
            if item not in itemCount.keys():
                itemCount[item] = 0

        pprint(itemCount)
        print items

        #now, populate combo box.
        if self.initDict == None:
            for item in items:
                itemLongName = self.dpdict[item]['longName']
                #self.qualifierCombo.addItem(itemLongName,item)

                rule = self.dpdict[item]['optional']
                if rule == 'one':
                    print 'one'
                    #Must have exactly one of this type of qualifier
                    if itemCount[item] == 0:
                        #add it to combo, mark as required.
                        self.dataPointCombo.addItem(itemLongName+' (REQUIRED)',item)
                    elif itemCount[item] == 0:
                        pass
                    else:
                        print 'More than one '+item+' - exactly one required.'
                        print 'XML is screwy - fix manually.'
                        exit()

                elif rule == 'zeroOrMore':
                    print 'zeroOrMore'
                    #can have as many or as few as we want, including none.
                    self.dataPointCombo.addItem(itemLongName,item)
                elif rule == 'oneOrMore':
                    print 'oneOrMore'
                    if itemCount[item] == 0:
                        #add it to combo, mark as required.
                        self.dataPointCombo.addItem(itemLongName+' (REQUIRED)',item)
                    elif itemCount[item] > 0:
                        self.dataPointCombo.addItem(itemLongName,item)
                    else:
                        print 'Something wrong with itemCount['+item+'] - should be'
                        print 'one or more.'
                        exit()

                elif rule == 'optional':
                    print 'optional'
                    if itemCount[item] == 0:
                        self.dataPointCombo.addItem(itemLongName,item)
                    elif itemCount[item] == 1:
                        pass
                    else:
                        print 'Something wrong with itemCount['+item+'] - should be'
                        print 'zero or one.'
                        exit()
                else:
                    print 'optional entry in dpdict is screwy. (CHEMICAL)',rule
                    exit()

        
        else:
            #editing a datapoint.
            item = self.initDict['DP_L0_DOMAIN_CODE']
            itemLongName = self.dpdict[item]['longName']
            self.dataPointCombo.addItem(itemLongName,item)



        #self.dataPointCombo.addItems(items)
        #sets index to current index - initialises combo box
        #self.dataPointCombo.setCurrentIndex(items.index('CONCLUSION'))

        #parent vbox
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addWidget(self.dataPointCombo)

        #now, add a child vbox to contain the dynamic stuff.
        self.dynamicvbox = QtGui.QVBoxLayout()
        self.vbox.addLayout(self.dynamicvbox)
        self.connect(self.dataPointCombo,QtCore.SIGNAL('currentIndexChanged(int)'),self.selectDataType)
        
        #make hbox for ok and cancel buttons.
        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        self.okBtn = QtGui.QPushButton('OK')
        self.okBtn.clicked.connect(self.accept)
        self.cancelBtn = QtGui.QPushButton('Cancel')
        self.cancelBtn.clicked.connect(self.reject)
        hbox.addWidget(self.okBtn)
        hbox.addWidget(self.cancelBtn)
        self.vbox.addLayout(hbox)

        self.setLayout(self.vbox)
        #self.setGeometry(50,50,300,300)
        self.dataPointCombo.setFocus()

        if self.initDict != None:
            oldIdx = self.dataPointCombo.currentIndex()
            self.dataPointCombo.setCurrentIndex(self.dataPointCombo.findData(QtCore.QString(self.initDict['DP_L0_DOMAIN_CODE']),QtCore.Qt.UserRole))
            if self.dataPointCombo.currentIndex() == oldIdx:
                #call selectDataPointType explicitly
                self.selectDataType(0)
        else:
            self.selectDataType(0)

    def selectDataType(self,selectedIdx):
        """Modifies dialog depending on what data type is chosen in combo box."""

        print 'selectDataType'
        if hasattr(self,'dynamicvbox'):    
            self.delLayout(self.dynamicvbox)

        self.dp = str(self.dataPointCombo.itemData(selectedIdx,QtCore.Qt.UserRole).toPyObject())
        self.dt = self.dpdict[str(self.dp)]['dpType']
        print 'Data Type (self.dpdict[str(self.dp)][\'dpType\']): ',[str(self.dp),self.dt]
        print self.dpdict[str(self.dp)]

        #print self.dt

        #need a dict of possible data_point names, their
        #types and units (if relevant).

        #now, add stuff pertinent to particular datatype
        if self.dt == 'DP_FLOAT':

            self.radioGroup = QtGui.QButtonGroup()

            self.singleValueRadioBtn = QtGui.QRadioButton('Single Value')
            self.rangeRadioBtn = QtGui.QRadioButton('Range')
            self.minRadioBtn = QtGui.QRadioButton('Min.')
            self.maxRadioBtn = QtGui.QRadioButton('Max.')

            if 'exceptions' in self.dpdict[str(self.dp)]:
                self.exceptRadioBtn = QtGui.QRadioButton('Exception')
            self.radioGroup.addButton(self.singleValueRadioBtn)
            self.radioGroup.addButton(self.rangeRadioBtn)
            self.radioGroup.addButton(self.minRadioBtn)
            self.radioGroup.addButton(self.maxRadioBtn)
            if 'exceptions' in self.dpdict[str(self.dp)]:
                self.radioGroup.addButton(self.exceptRadioBtn)
            self.radioGroup.buttonClicked.connect(self.floatTypeChanged)

            self.hbox1 = QtGui.QHBoxLayout()
            self.hbox1.addWidget(self.singleValueRadioBtn)
            self.hbox1.addWidget(self.rangeRadioBtn)
            self.hbox1.addWidget(self.minRadioBtn)
            self.hbox1.addWidget(self.maxRadioBtn)
            if 'exceptions' in self.dpdict[str(self.dp)]:
                self.hbox1.addWidget(self.exceptRadioBtn)
            self.dynamicvbox.addLayout(self.hbox1)

            #If self.initDict exists, work out how the float should
            #be displayed for editing (single value, range, min. value etc.)
            if self.initDict != None:
                if 'DP_QUANT_VALUE1' in self.initDict.keys() and \
                    'DP_QUANT_VALUE2' in self.initDict.keys():
                    #It's a range.
                    self.rangeRadioBtn.click()
                elif 'DP_QUANT_VALUE1' not in self.initDict.keys() and \
                    'DP_QUANT_VALUE2' in self.initDict.keys():
                    #Max value.
                    self.maxRadioBtn.click()
                elif 'DP_QUANT_VALUE1' in self.initDict.keys() and \
                    'DP_QUANT_VALUE2' not in self.initDict.keys():
                    #Single value or min value
                    if 'DP_QUAL_VALUE1' not in self.initDict.keys():
                        #single value
                        self.singleValueRadioBtn.click()
                    else:
                        #min value
                        self.minRadioBtn.click()
                elif 'DP_QUAL_VALUE1' in  self.initDict.keys() \
                        and 'DP_QUANT_VALUE1' not in self.initDict.keys()\
                        and 'DP_QUANT_VALUE2' not in self.initDict.keys()\
                        and 'exceptions' in self.dpdict[str(self.dp)]:
                    #exception
                    self.exceptRadioBtn.click()

                else:
                    print 'dataDialog: Don\'t know input data type.'
                    exit()
            else:
                self.singleValueRadioBtn.click()



        elif self.dt == 'DP_VARCHAR2':
            #print 'CLOB!!!!!',self.initDict['DP_CLOB1']
            self.rangeFormLayout = QtGui.QFormLayout()
            self.DP_QUAL_VALUE1 = QtGui.QTextEdit()
            #self.DP_QUAL_VALUE1.setFontPointSize(14)
            self.rangeFormLayout.addRow('Text input:',self.DP_QUAL_VALUE1)
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            #self.DP_QUAL_VALUE1.setFocus()
            
            #Initialise values if self.initDict exists
            if self.initDict != None:
                self.DP_QUAL_VALUE1.setPlainText(self.initDict['DP_QUAL_VALUE1'])
                try: 
                    if self.initDict['DP_CONFIDENTIAL1_YN'] == 'Y':
                        self.DP_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['DP_CONFIDENTIAL1_YN'] == 'N':
                        self.DP_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'DP_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass

        elif self.dt == 'DP_CLOB':
            #print 'CLOB!!!!!',self.initDict['DP_CLOB1']
            self.rangeFormLayout = QtGui.QFormLayout()
            self.DP_CLOB1 = QtGui.QTextEdit()
            #self.DP_CLOB1.setFontPointSize(14)
            #print self.DP_CLOB1.currentFont().key()
            self.rangeFormLayout.addRow('Text input:',self.DP_CLOB1)
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            #self.DP_CLOB1.setFocus()
            
            #Initialise values if self.initDict exists
            if self.initDict != None:
                self.DP_CLOB1.setPlainText(self.initDict['DP_CLOB1'])
                try: 
                    if self.initDict['DP_CONFIDENTIAL1_YN'] == 'Y':
                        self.DP_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['DP_CONFIDENTIAL1_YN'] == 'N':
                        self.DP_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'DP_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass

        elif self.dt == 'DP_YN':
            self.radioGroup =QtGui.QButtonGroup()
            self.yesRadioBtn = QtGui.QRadioButton('Yes')
            self.noRadioBtn = QtGui.QRadioButton('No')
            self.uncertainRadioBtn = QtGui.QRadioButton('Uncertain')

            self.radioGroup.addButton(self.yesRadioBtn)
            self.radioGroup.addButton(self.noRadioBtn)
            self.radioGroup.addButton(self.uncertainRadioBtn)

            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.yesRadioBtn)
            vbox1.addWidget(self.noRadioBtn)
            vbox1.addWidget(self.uncertainRadioBtn)
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            #self.yesRadioBtn.setFocus()

            #Initialise values if self.initDict exists
            if self.initDict != None:
                try:
                    if self.initDict['DP_QUAL_VALUE1'] == 'Y':
                        self.yesRadioBtn.setChecked(True)
                    elif self.initDict['DP_QUAL_VALUE1'] == 'N':
                        self.noRadioBtn.setChecked(True)
                    elif self.initDict['DP_QUAL_VALUE1'] == 'Uncertain':
                        self.uncertainRadioBtn.setChecked(True)
                    else:
                        print 'self.initDict[\'DP_QUAL_VALUE1\'] needs to be Yes or No or Uncertain.'
                except:
                    pass
                try:
                    if self.initDict['DP_CONFIDENTIAL1_YN'] == 'Y':
                        self.DP_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['DP_CONFIDENTIAL1_YN'] == 'N':
                        self.DP_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'DP_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass
            else:
                #make a default initial selection.
                self.yesRadioBtn.setChecked(True)
                
            
        
        elif self.dt == 'DP_VARCHAR2CHOICE':
            self.DP_QUAL_VALUE1 = QtGui.QComboBox()
            self.DP_QUAL_VALUE1.addItems(list(\
                    self.dpdict[str(self.dp)]['choices']))


            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.DP_QUAL_VALUE1)
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            #self.DP_QUAL_VALUE1.setFocus()


        elif self.dt == 'DP_INTEGER':
            self.DP_INTEGER1 = QtGui.QLineEdit()
            v = QtGui.QIntValidator()
            self.DP_INTEGER1.setValidator(v)
            
            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.DP_INTEGER1)
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            #self.DP_INTEGER1.setFocus()
            
            if self.initDict != None:
                self.DP_INTEGER1.setText(self.initDict['DP_INTEGER1'])
                try:
                    if self.initDict['DP_CONFIDENTIAL1_YN'] == 'Y':
                        self.DP_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['DP_CONFIDENTIAL1_YN'] == 'N':
                        self.DP_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'DP_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass
                try:
                    self.DP_UNIT_CODE1.setText(self.initDict['DP_UNIT_CODE1'])
                except:
                    pass
        
        else:
            print 'data type '+self.dt+' not valid' 
    


    def floatTypeChanged(self):
        """Changes float subtype depending on radio button choice."""
        #get rid of everything except the first item in 
        #in the dynamic layout. This is the box of radio buttons.
        for index in range(self.dynamicvbox.count()-1,0,-1):
            self.delLayoutItem(self.dynamicvbox,index)
        
        #make the relevant input fields.
        if self.radioGroup.checkedButton() == self.rangeRadioBtn:
            print 'rangeRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()
            
            self.DP_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.DP_QUANT_VALUE1.setValidator(dv)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE1 = QtGui.QComboBox()
                self.DP_UNIT_CODE1.addItems(list(\
                        self.dpdict[str(self.dp)]['dpUnits']))
            self.DP_QUAL_VALUE1 = QtGui.QComboBox()
            self.DP_QUAL_VALUE1.addItems(['>','>='])

            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Min. Confidential')
            
            self.DP_QUANT_VALUE2 = QtGui.QLineEdit()
            self.DP_QUANT_VALUE2.setValidator(dv)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE2 = QtGui.QComboBox()
                self.DP_UNIT_CODE2.addItems(list(\
                        self.dpdict[str(self.dp)]['dpUnits']))
            self.DP_QUAL_VALUE2 = QtGui.QComboBox()
            self.DP_QUAL_VALUE2.addItems(['<','<='])
            
            self.DP_CONFIDENTIAL2_YN = QtGui.QCheckBox('Max. Confidential')

            self.rangeFormLayout.addRow('> / >=',self.DP_QUAL_VALUE1)
            self.rangeFormLayout.addRow('Min. value',self.DP_QUANT_VALUE1)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.rangeFormLayout.addRow('Min. value unit',self.DP_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.rangeFormLayout.addRow('< / <=',self.DP_QUAL_VALUE2)
            self.rangeFormLayout.addRow('Max. value',self.DP_QUANT_VALUE2)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.rangeFormLayout.addRow('Max. value unit',self.DP_UNIT_CODE2)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL2_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)


            #make sure that if one unit is changed, the other is 
            #changed as well.
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE1.currentIndexChanged.connect\
                        (self.onDPUnitCode1Changed)
                self.DP_UNIT_CODE2.currentIndexChanged.connect\
                        (self.onDPUnitCode2Changed)
            
            self.DP_QUANT_VALUE1.setFocus()


        elif self.radioGroup.checkedButton() == self.minRadioBtn:
            print 'minRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()

            self.DP_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.DP_QUANT_VALUE1.setValidator(dv)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE1 = QtGui.QComboBox()
                self.DP_UNIT_CODE1.addItems(list(\
                        self.dpdict[str(self.dp)]['dpUnits']))
            self.DP_QUAL_VALUE1 = QtGui.QComboBox()
            self.DP_QUAL_VALUE1.addItems(['>','>='])
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('> / >=',self.DP_QUAL_VALUE1)
            self.rangeFormLayout.addRow('Min. value',self.DP_QUANT_VALUE1)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.rangeFormLayout.addRow('Min. value unit',self.DP_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
 
            self.DP_QUANT_VALUE1.setFocus()


        elif self.radioGroup.checkedButton() == self.maxRadioBtn:
            print 'maxRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()
            self.DP_QUANT_VALUE2 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.DP_QUANT_VALUE2.setValidator(dv)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE2 = QtGui.QComboBox()
                self.DP_UNIT_CODE2.addItems(list(\
                        self.dpdict[str(self.dp)]['dpUnits']))
            self.DP_QUAL_VALUE2 = QtGui.QComboBox()
            self.DP_QUAL_VALUE2.addItems(['<','<='])
            self.DP_CONFIDENTIAL2_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('< / <=',self.DP_QUAL_VALUE2)
            self.rangeFormLayout.addRow('Max. value',self.DP_QUANT_VALUE2)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.rangeFormLayout.addRow('Max. value unit',self.DP_UNIT_CODE2)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL2_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
 
            self.DP_QUANT_VALUE2.setFocus()

        
        elif self.radioGroup.checkedButton() == self.singleValueRadioBtn:
            self.rangeFormLayout = QtGui.QFormLayout()
            self.DP_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.DP_QUANT_VALUE1.setValidator(dv)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.DP_UNIT_CODE1 = QtGui.QComboBox()
                self.DP_UNIT_CODE1.addItems(list(\
                        self.dpdict[str(self.dp)]['dpUnits']))
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('Value',self.DP_QUANT_VALUE1)
            if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                self.rangeFormLayout.addRow('Value unit',self.DP_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            

        elif self.radioGroup.checkedButton() == self.exceptRadioBtn:
            self.rangeFormLayout = QtGui.QFormLayout()
            self.DP_QUAL_VALUE1 = QtGui.QComboBox()
            self.DP_QUAL_VALUE1.addItems(self.dpdict[str(self.dp)]['exceptions'])
            self.DP_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.DP_QUAL_VALUE1)
            self.rangeFormLayout.addRow(self.DP_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)

            self.DP_QUAL_VALUE1.setFocus()
         
        #Initialise values if self.initDict exists
        if self.initDict != None:
            try:
                self.DP_QUANT_VALUE1.setText(self.initDict['DP_QUANT_VALUE1'])
            except:
                pass
            try:
                self.DP_QUANT_VALUE2.setText(self.initDict['DP_QUANT_VALUE2'])
            except:
                pass
            try:
                self.DP_QUAL_VALUE1.setText(self.initDict['DP_QUAL_VALUE1'])
            except:
                pass
            try:
                self.DP_QUAL_VALUE2.setText(self.initDict['DP_QUAL_VALUE2'])
            except:
                pass
            try:
                if self.initDict['DP_CONFIDENTIAL1_YN'] == 'Y':
                    self.DP_CONFIDENTIAL1_YN.setChecked(True)
                elif self.initDict['DP_CONFIDENTIAL1_YN'] == 'N':
                    self.DP_CONFIDENTIAL1_YN.setChecked(False)
                else:
                    print 'self.initDict[\'DP_CONFIDENTIAL1\'] needs to be Y or N.'
            except:
                 pass
            try:
                if self.initDict['DP_CONFIDENTIAL2'] == 'Y':
                    self.DP_CONFIDENTIAL2.setChecked(True)
                elif self.initDict['DP_CONFIDENTIAL2'] == 'N':
                    self.DP_CONFIDENTIAL2.setChecked(False)
                else:
                    print 'self.initDict[\'DP_CONFIDENTIAL2\'] needs to be Y or N.'
            except:
                pass
            try:
                self.DP_UNIT_CODE1.setText(self.initDict['DP_UNIT_CODE1'])
            except:
                pass
            try:
                self.DP_UNIT_CODE2.setText(self.initDict['DP_UNIT_CODE2'])
            except:
                pass
               
               

    def onDPUnitCode1Changed(self):
        """Slot to change DP_UNIT_CODE2 selection when DP_UNIT_CODE1 seletion changes."""
        DPUnitCode1CurrentText  = self.DP_UNIT_CODE1.currentText()
        idx=self.DP_UNIT_CODE2.findText(DPUnitCode1CurrentText,\
                QtCore.Qt.MatchExactly)
        self.DP_UNIT_CODE2.setCurrentIndex(idx)

    def onDPUnitCode2Changed(self):
        """Slot to change DP_UNIT_CODE1 selection when DP_UNIT_CODE2 seletion changes."""
        DPUnitCode2CurrentText  = self.DP_UNIT_CODE2.currentText()
        idx=self.DP_UNIT_CODE1.findText(DPUnitCode2CurrentText,\
                QtCore.Qt.MatchExactly)
        self.DP_UNIT_CODE1.setCurrentIndex(idx)

    
    def delLayout(self,lo):
        """Clears the dynamic layout."""
        for index in range(lo.count()-1,-1,-1):
            self.delLayoutItem(lo,index)
    
    def delLayoutItem(self,lo,index):
        """Clears item in dynamic layout."""
        if issubclass(type(lo.itemAt(index)),QtGui.QLayout):
            item = lo.takeAt(index)
            self.delLayout(item)
        else:
            #print 'item: ',lo.itemAt(index)
            item = lo.takeAt(index)
            item.widget().deleteLater()
    def getData(self):
        """Gets data from dialog (once dialog has been closed)."""
        #first, need to work out what sort of datapoint it is and the name.
        print 'DATA_TYPE: ',self.dt

        dataout = {}
        dataout['DP_L0_DOMAIN_CODE'] = self.dp
        dataout['DP_L1_KINGDOM_CODE'] = self.dt

        if self.dt == 'DP_FLOAT':
            #This is the most complicated one.
            #Do stuff based on radio button status
            if self.radioGroup.checkedButton() == self.singleValueRadioBtn:
                dataout['DP_QUANT_VALUE1'] = str(self.DP_QUANT_VALUE1.text())
                print 'self.dpdict[str(self.dp)]: '
                print self.dpdict[str(self.dp)]
                if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                    dataout['DP_UNIT_CODE1'] = str(self.DP_UNIT_CODE1.currentText())

                if self.DP_CONFIDENTIAL1_YN.isChecked():
                    dataout['DP_CONFIDENTIAL1_YN']='Y'
                else:
                    dataout['DP_CONFIDENTIAL1_YN']='N'

            if self.radioGroup.checkedButton() == self.rangeRadioBtn or self.radioGroup.checkedButton() == self.minRadioBtn :
                dataout['DP_QUAL_VALUE1'] = str(self.DP_QUAL_VALUE1.currentText())
                dataout['DP_QUANT_VALUE1'] = str(self.DP_QUANT_VALUE1.text())
                if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                    dataout['DP_UNIT_CODE1'] = str(self.DP_UNIT_CODE1.currentText())
                if self.DP_CONFIDENTIAL1_YN.isChecked():
                    dataout['DP_CONFIDENTIAL1_YN']='Y'
                else:
                    dataout['DP_CONFIDENTIAL1_YN']='N'

            if self.radioGroup.checkedButton() == self.rangeRadioBtn or self.radioGroup.checkedButton() == self.maxRadioBtn :
                dataout['DP_QUAL_VALUE2'] = str(self.DP_QUAL_VALUE2.currentText())
                dataout['DP_QUANT_VALUE2'] = str(self.DP_QUANT_VALUE2.text())
                if len(self.dpdict[str(self.dp)]['dpUnits']) > 0:
                    dataout['DP_UNIT_CODE2'] = str(self.DP_UNIT_CODE2.currentText())
                if self.DP_CONFIDENTIAL2_YN.isChecked():
                    dataout['DP_CONFIDENTIAL2_YN']='Y'
                else:
                    dataout['DP_CONFIDENTIAL2_YN']='N'
            if hasattr(self,'exceptRadioBtn'):
                if self.radioGroup.checkedButton() == self.exceptRadioBtn:
                    dataout['DP_QUAL_VALUE1'] = str(self.DP_QUAL_VALUE1.currentText())
                    if self.DP_CONFIDENTIAL1_YN.isChecked():
                        dataout['DP_CONFIDENTIAL1_YN']='Y'
                    else:
                        dataout['DP_CONFIDENTIAL1_YN']='N'

        elif self.dt == 'DP_VARCHAR2':
            dataout['DP_QUAL_VALUE1'] = str(self.DP_QUAL_VALUE1.toPlainText())
            if self.DP_CONFIDENTIAL1_YN.isChecked():
                dataout['DP_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['DP_CONFIDENTIAL1_YN']='N'

        elif self.dt == 'DP_CLOB':
            dataout['DP_CLOB1'] = str(self.DP_CLOB1.toPlainText())
            if self.DP_CONFIDENTIAL1_YN.isChecked():
                dataout['DP_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['DP_CONFIDENTIAL1_YN']='N'

        elif self.dt == 'DP_YN':
            if self.radioGroup.checkedButton() == self.yesRadioBtn:
                dataout['DP_QUAL_VALUE1'] = 'Y'
            elif self.radioGroup.checkedButton() == self.noRadioBtn:
                dataout['DP_QUAL_VALUE1'] = 'N'
            elif self.radioGroup.checkedButton() == self.uncertainRadioBtn:
                dataout['DP_QUAL_VALUE1'] = 'Uncertain'
            else:
                print 'DP_YN: ???'

            if self.DP_CONFIDENTIAL1_YN.isChecked():
                dataout['DP_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['DP_CONFIDENTIAL1_YN']='N'

        elif self.dt == 'DP_VARCHAR2CHOICE':
            if str(self.DP_QUAL_VALUE1.currentText()) != '':
                dataout['DP_QUAL_VALUE1'] = str(self.DP_QUAL_VALUE1.currentText())
            else:
                print 'self.DP_QUAL_VALUE1 has value \'\''
                exit()
            if self.DP_CONFIDENTIAL1_YN.isChecked():
                dataout['DP_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['DP_CONFIDENTIAL1_YN']='N'

        elif self.dt == 'DP_INTEGER':
            print 'self.DP_INTEGER1.text():' + str(self.DP_INTEGER1.text())
            dataout['DP_INTEGER1'] = str(self.DP_INTEGER1.text())
            
            if self.DP_CONFIDENTIAL1_YN.isChecked():
                dataout['DP_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['DP_CONFIDENTIAL1_YN']='N'

        return dataout

    def accept(self):
        """Re-implements accept to check for empty fields."""
        dataout = self.getData()
        emptyFlag = False

        if dataout['DP_L1_KINGDOM_CODE'] == 'DP_FLOAT':
            #check to see whether there are blank input fields
            print 'It\'s a float!'
            try:
                print 'DP_QUANT_VALUE1: ',[dataout['DP_QUANT_VALUE1']]
            except:
                pass
            try:
                print 'DP_QUANT_VALUE2',[dataout['DP_QUANT_VALUE2']]
            except:
                pass

            if 'DP_QUANT_VALUE1' in dataout.keys():
                if dataout['DP_QUANT_VALUE1'] == '':
                    emptyFlag = True
            if 'DP_QUANT_VALUE2' in dataout.keys():
                if dataout['DP_QUANT_VALUE2'] == '':
                    emptyFlag = True

        elif dataout['DP_L1_KINGDOM_CODE'] == 'DP_CLOB':
            if 'DP_CLOB1' in dataout.keys():
                if dataout['DP_CLOB1'].strip() == '':
                    emptyFlag = True
        elif dataout['DP_L1_KINGDOM_CODE'] == 'DP_VARCHAR2':
            if 'DP_QUAL_VALUE1' in dataout.keys():
                if dataout['DP_QUAL_VALUE1'].strip() == '':
                    emptyFlag = True

        if emptyFlag == True:
            QtGui.QMessageBox.information(self,'Cascading Viewer',\
                'Can\'t have empty fields!', QtGui.QMessageBox.Ok)
        else:
            return super(dataDialog,self).accept()
                    

class qualifierDialog(QtGui.QDialog):
    """Dialog box to allow addition of a new qualifier.
    
    This is also used to make a copy of an old qualifier
    for further editing."""
    def __init__(self,parent,dpDict,initDict=None):
        super(qualifierDialog,self).__init__(parent)
        if initDict !=None:
            self.initDict = initDict
        else:
            self.initDict = None

        self.dpDict = dpDict

        print 'qualifierDialog.initDict: ', self.initDict

        #self.schemaXML = etree.parse('newschema.rng').getroot()
        #self.schemaRNG = etree.RelaxNG(self.schemaXML)

        #self.dpDict = self.dictFromSchema()

        #need to know about the parent to be able to
        #work out what qualifiers are allowed for that
        #particular data point
        #print 'parent: ',dir(parent)
        self.parentDpIdx = parent.dataPointProxyModel.mapToSource(parent.dataPointSelectionModel.currentIndex())
        self.parentTag = str(self.parentDpIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'DP_L0_DOMAIN_CODE')].toString())

        #print self.parentDpIdx


        self.initUI()

    def initUI(self):
        """Initialises user interface of dialog."""
       
        #make combo box.
        self.qualifierCombo = QtGui.QComboBox()
        #self.dpdict = self.dictFromSchema()

        #list of all possible entries
        items = self.dpDict[self.parentTag]['qualifiers'].keys()
        
        #count all of the items the datapoint currently has.
        itemCount = {}
        for item in self.parentDpIdx.internalPointer().childItems:
            itemName = item.datadict['DP_L0_DOMAIN_CODE']
            if itemName not in itemCount.keys():
                itemCount[itemName] = 1
            else:
                itemCount[itemName] += 1

        #For everything else in items, set itemCount to zero.
        for item in items:
            if item not in itemCount.keys():
                itemCount[item] = 0
        
        pprint(itemCount)
        print items

        #Now, populate combo box.
        if self.initDict == None:
            for item in items:
                itemLongName = self.dpDict[self.parentTag]['qualifiers'][item]['longName']
                #self.qualifierCombo.addItem(itemLongName,item)

                rule = self.dpDict[self.parentTag]['qualifiers'][item]['optional']
                if rule == 'one':
                    print 'one'
                    #Must have exactly one of this type of qualifier
                    if itemCount[item] == 0:
                        #add it to combo, mark as required.
                        self.qualifierCombo.addItem(itemLongName+' (REQUIRED)',item)
                    elif itemCount[item] == 0:
                        pass
                    else:
                        print 'More than one '+item+' - exactly one required.'
                        print 'XML is screwy - fix manually.'
                        exit()

                elif rule == 'zeroOrMore':
                    print 'zeroOrMore'
                    #can have as many or as few as we want, including none.
                    self.qualifierCombo.addItem(itemLongName,item)
                elif rule == 'oneOrMore':
                    print 'oneOrMore'
                    if itemCount[item] == 0:
                        #add it to combo, mark as required.
                        self.qualifierCombo.addItem(itemLongName+' (REQUIRED)',item)
                    elif itemCount[item] > 0:
                        self.qualifierCombo.addItem(itemLongName,item)
                    else:
                        print 'Something wrong with itemCount['+item+'] - should be'
                        print 'one or more.'
                        exit()
                   
                elif rule == 'optional':
                    print 'optional'
                    if itemCount[item] == 0:
                        self.qualifierCombo.addItem(itemLongName,item)
                    elif itemCount[item] == 1:
                        pass
                    else:
                        print 'Something wrong with itemCount['+item+'] - should be'
                        print 'zero or one.'
                        exit()

                else:
                    print 'optional entry in dpDict is screwy.',rule,self.parentTag,qName
                    exit()
        else:
            #editing a qualifier. #get input name
            print self.initDict
            item = self.initDict['DP_L0_DOMAIN_CODE']
            itemLongName = self.dpDict[self.parentTag]['qualifiers'][item]['longName']
            self.qualifierCombo.addItem(itemLongName,item)

        #print items
        #self.qualifierCombo.addItems(items)
        #sets index to current index - initialises combo box
        #self.dataPointCombo.setCurrentIndex(items.index('CONCLUSION'))

        #parent vbox
        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addWidget(self.qualifierCombo)

        #now, add a child vbox to contain the dynamic stuff.
        self.dynamicvbox = QtGui.QVBoxLayout()
        self.vbox.addLayout(self.dynamicvbox)
        self.connect(self.qualifierCombo,QtCore.SIGNAL('currentIndexChanged(int)'),self.selectQualifierType)

        #Note that if the qualifier you are editing is the first one
        #in the combo, and you set the current index, then currentIndexCHanged is
        #not emitted. In this case, call selectQualifierType explicitly
        #if self.initDict != None:
        #    oldIdx = self.qualifierCombo.currentIndex()

        #    print 'DP_L0_DOMAIN_CODE: ',self.initDict['DP_L0_DOMAIN_CODE']
        #    self.qualifierCombo.setCurrentIndex(self.qualifierCombo.findData(QtCore.QString(self.initDict['DP_L0_DOMAIN_CODE']),QtCore.Qt.DisplayRole))
        #    if self.qualifierCombo.currentIndex() == oldIdx:
        #        #call selectQualifierType explicitly
        #        self.selectQualifierType(self.initDict['DP_L0_DOMAIN_CODE'])
        #else:
        #    self.qualifierCombo.setCurrentIndex(1)
        #    self.qualifierCombo.setCurrentIndex(0)


        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        self.okBtn = QtGui.QPushButton('OK')
        self.okBtn.clicked.connect(self.accept)
        self.cancelBtn = QtGui.QPushButton('Cancel')
        self.cancelBtn.clicked.connect(self.reject)
        hbox.addWidget(self.okBtn)
        hbox.addWidget(self.cancelBtn)
        self.vbox.addLayout(hbox)

        self.setLayout(self.vbox)
        #self.setGeometry(50,50,300,300)
        self.qualifierCombo.setFocus()
        print 'done qualifierDialog.initUI'
        
        if self.initDict != None:
            oldIdx = self.qualifierCombo.currentIndex()
            self.qualifierCombo.setCurrentIndex(self.qualifierCombo.findData(QtCore.QString(self.initDict['DP_L0_DOMAIN_CODE']),QtCore.Qt.UserRole))

            if self.qualifierCombo.currentIndex() == oldIdx:
                #call selectQualifierType explicitly
                self.selectQualifierType(oldIdx)
        else:
            #New qualifier.
            #self.selectQualifierType(str(self.qualifierCombo.itemData(0,QtCore.Qt.UserRole).toPyObject()))
            self.selectQualifierType(0)
 
    def selectQualifierType(self,selectedIdx):
        """Modifies the dialog based on the data type of the selected qualifier."""
        print 'selectQualifierType',selectedIdx

        if hasattr(self,'dynamicvbox'):
            self.delLayout(self.dynamicvbox)

        self.rangeFormLayout = QtGui.QFormLayout()

        self.qName = str(self.qualifierCombo.itemData(selectedIdx,QtCore.Qt.UserRole).toPyObject())
        print 'self.qName: ',self.qName
        #self.qName = str(qName)
        self.qt = self.dpDict[self.parentTag]['qualifiers'][str(self.qName)]['qType']
        print 'Qualifier name, Qualifier type: ',str(self.qName),self.qt
        pprint(self.dpDict[self.parentTag]['qualifiers'])
        

        if self.qt == 'Q_TAG':
            print 'It\'s a Q_TAG!'
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)

        elif self.qt == 'Q_CLOB':
            print 'It\'s a Q_CLOB!'

            self.Q_CLOB1 = QtGui.QTextEdit()
            #self.Q_CLOB1.setFontPointSize(14)
            self.rangeFormLayout.addRow('Text input:',self.Q_CLOB1)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
           # self.Q_CLOB1.setFocus()
            
            #Initialise values if self.initDict exists
            if self.initDict != None:
                self.Q_CLOB1.setHtml(self.initDict['Q_CLOB1'])
                try: 
                    if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                        self.Q_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                        self.Q_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass

        elif self.qt == 'Q_VARCHAR2':
            print 'It\'s a Q_VARCHAR2!'

            self.Q_QUAL_VALUE1 = QtGui.QTextEdit()
            #self.Q_QUAL_VALUE1.setFontPointSize(14)
            self.rangeFormLayout.addRow('Text input:',self.Q_QUAL_VALUE1)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            #self.Q_QUAL_VALUE1.setFocus()
            
            #Initialise values if self.initDict exists
            if self.initDict != None:
                self.Q_QUAL_VALUE1.setHtml(self.initDict['Q_QUAL_VALUE1'])
                try: 
                    if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                        self.Q_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                        self.Q_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass

        elif self.qt == 'Q_FLOAT':

            self.radioGroup = QtGui.QButtonGroup()

            self.singleValueRadioBtn = QtGui.QRadioButton('Single Value')
            self.rangeRadioBtn = QtGui.QRadioButton('Range')
            self.minRadioBtn = QtGui.QRadioButton('Min.')
            self.maxRadioBtn = QtGui.QRadioButton('Max.')
            self.radioGroup.addButton(self.singleValueRadioBtn)
            self.radioGroup.addButton(self.rangeRadioBtn)
            self.radioGroup.addButton(self.minRadioBtn)
            self.radioGroup.addButton(self.maxRadioBtn)
            self.radioGroup.buttonClicked.connect(self.floatTypeChanged)

            self.hbox1 = QtGui.QHBoxLayout()
            self.hbox1.addWidget(self.singleValueRadioBtn)
            self.hbox1.addWidget(self.rangeRadioBtn)
            self.hbox1.addWidget(self.minRadioBtn)
            self.hbox1.addWidget(self.maxRadioBtn)
            self.dynamicvbox.addLayout(self.hbox1)
            self.dynamicvbox.addLayout(self.rangeFormLayout)

            #If self.initDict exists, work out how the float should
            #be displayed for editing (single value, range, min. value etc.)
            if self.initDict != None:
                if 'Q_QUANT_VALUE1' in self.initDict.keys() and \
                    'Q_QUANT_VALUE2' in self.initDict.keys():
                    #It's a range.
                    self.rangeRadioBtn.click()
                elif 'Q_QUANT_VALUE1' not in self.initDict.keys() and \
                    'Q_QUANT_VALUE2' in self.initDict.keys():
                    #Max value.
                    self.maxRadioBtn.click()
                elif 'Q_QUANT_VALUE1' in self.initDict.keys() and \
                    'Q_QUANT_VALUE2' not in self.initDict.keys():
                    #Single value or min value
                    if 'Q_QUAL_VALUE1' not in self.initDict.keys():
                        #single value
                        self.singleValueRadioBtn.click()
                    else:
                        #min value
                        self.minRadioBtn.click()
            else:
                self.singleValueRadioBtn.click()

        elif self.qt == 'Q_YN':
            self.radioGroup =QtGui.QButtonGroup()
            self.yesRadioBtn = QtGui.QRadioButton('Yes')
            self.noRadioBtn = QtGui.QRadioButton('No')

            self.radioGroup.addButton(self.yesRadioBtn)
            self.radioGroup.addButton(self.noRadioBtn)

            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.yesRadioBtn)
            vbox1.addWidget(self.noRadioBtn)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            self.dynamicvbox.addLayout(self.rangeFormLayout)

            #Initialise values if self.initDict exists
            if self.initDict != None:
                try:
                    if self.initDict['Q_QUAL_VALUE1'] == 'Y':
                        self.yesRadioBtn.setChecked(True)
                    elif self.initDict['Q_QUAL_VALUE1'] == 'N':
                        self.noRadioBtn.setChecked(True)
                    else:
                        print 'self.initDict[\'Q_QUAL_VALUE1\'] needs to be Y or N.'
                except:
                    pass
                try:
                    if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                        self.Q_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                        self.Q_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass
            else:
                #make a default initial selection.
                self.yesRadioBtn.setChecked(True)
            
        
        elif self.qt == 'Q_YNPStar':
            self.radioGroup =QtGui.QButtonGroup()
            self.yesRadioBtn = QtGui.QRadioButton('Yes')
            self.noRadioBtn = QtGui.QRadioButton('No')
            self.PStarRadioBtn = QtGui.QRadioButton('P*')

            self.radioGroup.addButton(self.yesRadioBtn)
            self.radioGroup.addButton(self.noRadioBtn)
            self.radioGroup.addButton(self.PStarRadioBtn)

            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.yesRadioBtn)
            vbox1.addWidget(self.noRadioBtn)
            vbox1.addWidget(self.PStarRadioBtn)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            self.dynamicvbox.addLayout(self.rangeFormLayout)

            #Initialise values if self.initDict exists
            if self.initDict != None:
                try:
                    if self.initDict['Q_QUAL_VALUE1'] == 'Y':
                        self.yesRadioBtn.setChecked(True)
                    elif self.initDict['Q_QUAL_VALUE1'] == 'N':
                        self.noRadioBtn.setChecked(True)
                    elif self.initDict['Q_QUAL_VALUE1'] == 'P*':
                        self.PStarRadioBtn.setChecked(True)
                    else:
                        print 'self.initDict[\'Q_QUAL_VALUE1\'] needs to be Y or N or P*.'
                except:
                    pass
                try:
                    if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                        self.Q_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                        self.Q_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass
                
            else:
                #make a default initial selection.
                self.yesRadioBtn.setChecked(True)
                
           
            

        elif self.qt == 'Q_INTEGER':
            self.Q_INTEGER1 = QtGui.QLineEdit()
            v = QtGui.QIntValidator()
            self.Q_INTEGER1.setValidator(v)
            
            vbox1 = QtGui.QVBoxLayout()
            vbox1.addWidget(self.Q_INTEGER1)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')
            vbox1.addWidget(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(vbox1)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            #self.Q_INTEGER1.setFocus()
            
            if self.initDict != None:
                self.Q_INTEGER1.setText(self.initDict['Q_INTEGER1'])
                try:
                    if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                        self.Q_CONFIDENTIAL1_YN.setChecked(True)
                    elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                        self.Q_CONFIDENTIAL1_YN.setChecked(False)
                    else:
                        print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
                except:
                    pass
                try:
                    self.Q_UNIT_CODE1.setText(self.initDict['Q_UNIT_CODE1'])
                except:
                    pass
               
   
        else:
            print 'self.qt '+self.qt+ 'is not a valid type.'
        
        #self.dynamicvbox.addLayout(self.rangeFormLayout)



    def floatTypeChanged(self):
        """Changes float subtype depending on radio button choice."""
        #get rid of everything except the first item in 
        #in the dynamic layout. This is the box of radio buttons.
        for index in range(self.dynamicvbox.count()-1,0,-1):
            self.delLayoutItem(self.dynamicvbox,index)


        #print 'Units: ',self.dpDict[self.parentTag]['qualifiers'][self.qName]
        
        #make the relevant input fields.
        if self.radioGroup.checkedButton() == self.rangeRadioBtn:
            print 'rangeRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()
            
            self.Q_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.Q_QUANT_VALUE1.setValidator(dv)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.Q_UNIT_CODE1 = QtGui.QComboBox()
                self.Q_UNIT_CODE1.addItems(list(\
                        self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']))
            self.Q_QUAL_VALUE1 = QtGui.QComboBox()
            self.Q_QUAL_VALUE1.addItems(['>','>='])

            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Min. Confidential')
            
            self.Q_QUANT_VALUE2 = QtGui.QLineEdit()
            self.Q_QUANT_VALUE2.setValidator(dv)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.Q_UNIT_CODE2 = QtGui.QComboBox()
                self.Q_UNIT_CODE2.addItems(list(\
                        self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']))
                        
            self.Q_QUAL_VALUE2 = QtGui.QComboBox()
            self.Q_QUAL_VALUE2.addItems(['<','<='])
            
            self.Q_CONFIDENTIAL2_YN = QtGui.QCheckBox('Max. Confidential')

            self.rangeFormLayout.addRow('> / >=',self.Q_QUAL_VALUE1)
            self.rangeFormLayout.addRow('Min. value',self.Q_QUANT_VALUE1)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.rangeFormLayout.addRow('Min. value unit',self.Q_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.rangeFormLayout.addRow('< / <=',self.Q_QUAL_VALUE2)
            self.rangeFormLayout.addRow('Max. value',self.Q_QUANT_VALUE2)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.rangeFormLayout.addRow('Max. value unit',self.Q_UNIT_CODE2)
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL2_YN)

            #make sure that if one unit is changed, the other is 
            #changed as well.
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.Q_UNIT_CODE1.currentIndexChanged.connect\
                        (self.onQUnitCode1Changed)
                self.Q_UNIT_CODE2.currentIndexChanged.connect\
                        (self.onQUnitCode2Changed)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            self.Q_QUANT_VALUE1.setFocus()
            


        elif self.radioGroup.checkedButton() == self.minRadioBtn:
            print 'minRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()

            self.Q_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.Q_QUANT_VALUE1.setValidator(dv)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.Q_UNIT_CODE1 = QtGui.QComboBox()
                self.Q_UNIT_CODE1.addItems(list(\
                        self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']))
            self.Q_QUAL_VALUE1 = QtGui.QComboBox()
            self.Q_QUAL_VALUE1.addItems(['>','>='])
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('> / >=',self.Q_QUAL_VALUE1)
            self.rangeFormLayout.addRow('Min. value',self.Q_QUANT_VALUE1)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.rangeFormLayout.addRow('Min. value unit',self.Q_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            self.Q_QUANT_VALUE1.setFocus()
 

        elif self.radioGroup.checkedButton() == self.maxRadioBtn:
            print 'maxRadioBtn'
            self.rangeFormLayout = QtGui.QFormLayout()
            self.Q_QUANT_VALUE2 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.Q_QUANT_VALUE2.setValidator(dv)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.Q_UNIT_CODE2 = QtGui.QComboBox()
                self.Q_UNIT_CODE2.addItems(list(\
                       self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits'] ))
            self.Q_QUAL_VALUE2 = QtGui.QComboBox()
            self.Q_QUAL_VALUE2.addItems(['<','<='])
            self.Q_CONFIDENTIAL2_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('< / <=',self.Q_QUAL_VALUE2)
            self.rangeFormLayout.addRow('Max. value',self.Q_QUANT_VALUE2)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.rangeFormLayout.addRow('Max. value unit',self.Q_UNIT_CODE2)
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL2_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            self.Q_QUANT_VALUE2.setFocus()
 
        
        elif self.radioGroup.checkedButton() == self.singleValueRadioBtn:
            self.rangeFormLayout = QtGui.QFormLayout()
            self.Q_QUANT_VALUE1 = QtGui.QLineEdit()
            #dv = QtGui.QDoubleValidator()
            dv = newDoubleValidator()
            self.Q_QUANT_VALUE1.setValidator(dv)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                tempUnits = list(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits'])
                if tempUnits == ['day','hour']:
                    tempUnits = ['hour','day']
                self.Q_UNIT_CODE1 = QtGui.QComboBox()
                self.Q_UNIT_CODE1.addItems(tempUnits)
            self.Q_CONFIDENTIAL1_YN = QtGui.QCheckBox('Confidential')

            self.rangeFormLayout.addRow('Value',self.Q_QUANT_VALUE1)
            if len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits']) > 0:
                self.rangeFormLayout.addRow('Value unit',self.Q_UNIT_CODE1)
            self.rangeFormLayout.addRow(self.Q_CONFIDENTIAL1_YN)
            self.dynamicvbox.addLayout(self.rangeFormLayout)
            #self.Q_QUANT_VALUE1.setFocus()

 
        #Initialise values if self.initDict exists
        if self.initDict != None:
            try:
                self.Q_QUANT_VALUE1.setText(self.initDict['Q_QUANT_VALUE1'])
            except:
                pass
            try:
                self.Q_QUANT_VALUE2.setText(self.initDict['Q_QUANT_VALUE2'])
            except:
                pass
            try:
                self.Q_QUAL_VALUE1.setText(self.initDict['Q_QUAL_VALUE1'])
            except:
                pass
            try:
                self.Q_QUAL_VALUE2.setText(self.initDict['Q_QUAL_VALUE2'])
            except:
                pass
            try:
                if self.initDict['Q_CONFIDENTIAL1_YN'] == 'Y':
                    self.Q_CONFIDENTIAL1_YN.setChecked(True)
                elif self.initDict['Q_CONFIDENTIAL1_YN'] == 'N':
                    self.Q_CONFIDENTIAL1_YN.setChecked(False)
                else:
                    print 'self.initDict[\'Q_CONFIDENTIAL1\'] needs to be Y or N.'
            except:
                pass
            try:
                if self.initDict['Q_CONFIDENTIAL2'] == 'Y':
                    self.Q_CONFIDENTIAL2.setChecked(True)
                elif self.initDict['Q_CONFIDENTIAL2'] == 'N':
                    self.Q_CONFIDENTIAL2.setChecked(False)
                else:
                    print 'self.initDict[\'Q_CONFIDENTIAL2\'] needs to be Y or N.'
            except:
                pass
            try:
                self.Q_UNIT_CODE1.setText(self.initDict['Q_UNIT_CODE1'])
            except:
                pass
            try:
                self.Q_UNIT_CODE2.setText(self.initDict['Q_UNIT_CODE2'])
            except:
                pass
               
        #self.dynamicvbox.addLayout(self.rangeFormLayout)
               

    def onQUnitCode1Changed(self):
        """Slot to change DP_UNIT_CODE2 selection when DP_UNIT_CODE1 seletion changes.
        
        The names of the combos are the same as in the data dialog 
        (where the code was pinched from)."""
        DPUnitCode1CurrentText  = self.DP_UNIT_CODE1.currentText()
        idx=self.DP_UNIT_CODE2.findText(DPUnitCode1CurrentText,\
                QtCore.Qt.MatchExactly)
        self.DP_UNIT_CODE2.setCurrentIndex(idx)

    def onQUnitCode2Changed(self):
        """Slot to change DP_UNIT_CODE1 selection when DP_UNIT_CODE2 seletion changes.
        
        The names of the combos are the same as in the data dialog 
        (where the code was pinched from)."""
        DPUnitCode2CurrentText  = self.DP_UNIT_CODE2.currentText()
        idx=self.DP_UNIT_CODE1.findText(DPUnitCode2CurrentText,\
                QtCore.Qt.MatchExactly)
        self.DP_UNIT_CODE1.setCurrentIndex(idx)



    def delLayout(self,lo):
        """Clears layout."""
        for index in range(lo.count()-1,-1,-1):
            self.delLayoutItem(lo,index)

    def delLayoutItem(self,lo,index):
        """Clears item in layout."""
        if issubclass(type(lo.itemAt(index)),QtGui.QLayout):
            item = lo.takeAt(index)
            self.delLayout(item)
        else:
            #print 'item: ',lo.itemAt(index)
            item = lo.takeAt(index)
            item.widget().deleteLater()

    def getData(self):
        """Gets data from the dialog once it has been closed."""
        dataout = {}
        #Note that both data_points and qualifiers have 
        #DP_L0_DOMAIN_CODE in their dicts - probably
        #should be Q_L0_DOMAIN_CODE.
        dataout['DP_L0_DOMAIN_CODE'] = self.qName
        dataout['Q_L1_KINGDOM_CODE'] = self.qt

        if self.qt == 'Q_TAG':
            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        if self.qt == 'Q_FLOAT':
            #This is the most complicated one.
            #Do stuff based on radio button status
            qUnitLen = len(self.dpDict[self.parentTag]['qualifiers'][self.qName]['qUnits'])
            if self.radioGroup.checkedButton() == self.singleValueRadioBtn:
                dataout['Q_QUANT_VALUE1'] = str(self.Q_QUANT_VALUE1.text())
                if qUnitLen > 0:
                    dataout['Q_UNIT_CODE1'] = str(self.Q_UNIT_CODE1.currentText())

                if self.Q_CONFIDENTIAL1_YN.isChecked():
                    dataout['Q_CONFIDENTIAL1_YN']='Y'
                else:
                    dataout['Q_CONFIDENTIAL1_YN']='N'

            if self.radioGroup.checkedButton() == self.rangeRadioBtn or self.radioGroup.checkedButton() == self.minRadioBtn :
                dataout['Q_QUAL_VALUE1'] = str(self.Q_QUAL_VALUE1.currentText())
                dataout['Q_QUANT_VALUE1'] = str(self.Q_QUANT_VALUE1.text())
                if qUnitLen > 0:
                    dataout['Q_UNIT_CODE1'] = str(self.Q_UNIT_CODE1.currentText())
                if self.Q_CONFIDENTIAL1_YN.isChecked():
                    dataout['Q_CONFIDENTIAL1_YN']='Y'
                else:
                    dataout['Q_CONFIDENTIAL1_YN']='N'

            if self.radioGroup.checkedButton() == self.rangeRadioBtn or self.radioGroup.checkedButton() == self.maxRadioBtn :
                dataout['Q_QUAL_VALUE2'] = str(self.Q_QUAL_VALUE2.currentText())
                dataout['Q_QUANT_VALUE2'] = str(self.Q_QUANT_VALUE2.text())
                if qUnitLen > 0:
                    dataout['Q_UNIT_CODE2'] = str(self.Q_UNIT_CODE2.currentText())
                if self.Q_CONFIDENTIAL2_YN.isChecked():
                    dataout['Q_CONFIDENTIAL2_YN']='Y'
                else:
                    dataout['Q_CONFIDENTIAL2_YN']='N'
        elif self.qt == 'Q_CLOB':
            dataout['Q_CLOB1'] = str(self.Q_CLOB1.toPlainText())
            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        elif self.qt == 'Q_VARCHAR2':
            dataout['Q_QUAL_VALUE1'] = str(self.Q_QUAL_VALUE1.toPlainText())
            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        elif self.qt == 'Q_YN':
            if self.radioGroup.checkedButton() == self.yesRadioBtn:
                dataout['Q_QUAL_VALUE1'] = 'Y'
            elif self.radioGroup.checkedButton() == self.noRadioBtn:
                dataout['Q_QUAL_VALUE1'] = 'N'
            else:
                print 'Q_YN???'

            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        elif self.qt == 'Q_YNPStar':
            if self.radioGroup.checkedButton() == self.yesRadioBtn:
                dataout['Q_QUAL_VALUE1'] = 'Y'
            elif self.radioGroup.checkedButton() == self.noRadioBtn:
                dataout['Q_QUAL_VALUE1'] = 'N'
            elif self.radioGroup.checkedButton() == self.PStarRadioBtn:
                dataout['Q_QUAL_VALUE1'] = 'P*'
            else:
                print 'Q_YNPStar???'
            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        elif self.qt == 'Q_INTEGER':
            print 'self.Q_INTEGER1.text():' + str(self.Q_INTEGER1.text())
            dataout['Q_INTEGER1'] = str(self.Q_INTEGER1.text())
            
            if self.Q_CONFIDENTIAL1_YN.isChecked():
                dataout['Q_CONFIDENTIAL1_YN']='Y'
            else:
                dataout['Q_CONFIDENTIAL1_YN']='N'

        return dataout

    def accept(self):
        """Re-implements accept to check for empty fields."""
        dataout = self.getData()
        emptyFlag = False

        if dataout['Q_L1_KINGDOM_CODE'] == 'Q_FLOAT':
            #check to see whether there are blank input fields
            print 'It\'s a float!'
            try:
                print 'Q_QUANT_VALUE1: ',[dataout['Q_QUANT_VALUE1']]
            except:
                pass
            try:
                print 'Q_QUANT_VALUE2',[dataout['Q_QUANT_VALUE2']]
            except:
                pass

            if 'Q_QUANT_VALUE1' in dataout.keys():
                if dataout['Q_QUANT_VALUE1'] == '':
                    emptyFlag = True
            if 'Q_QUANT_VALUE2' in dataout.keys():
                if dataout['Q_QUANT_VALUE2'] == '':
                    emptyFlag = True

        elif dataout['Q_L1_KINGDOM_CODE'] == 'Q_CLOB':
            if 'Q_CLOB1' in dataout.keys():
                if dataout['Q_CLOB1'].strip() == '':
                    emptyFlag = True
        elif dataout['Q_L1_KINGDOM_CODE'] == 'Q_VARCHAR2':
            if 'Q_QUAL_VALUE1' in dataout.keys():
                if dataout['Q_QUAL_VALUE1'].strip() == '':
                    emptyFlag = True

        if emptyFlag == True:
            QtGui.QMessageBox.information(self,'Cascading Viewer',\
                'Can\'t have empty fields!', QtGui.QMessageBox.Ok)
        else:
            return super(qualifierDialog,self).accept()
                   

class assessmentTool(QtGui.QMainWindow):
    """Contains cascadingViewer and provides the menu bar."""
    
    def __init__(self):
        
        super(assessmentTool,self).__init__()
        self.initData()
        self.initUI()

    def initData(self):
        """Initialises data."""

        self.schemaXML = etree.parse('newschema.rng').getroot()
        self.schemaRNG = etree.RelaxNG(self.schemaXML)

        #self.xml = etree.parse(self.xmlFileName)
        #self.xmlWorking = deepcopy(self.xml)
        ##xmlin = etree.ElementTree(etree.XML('<ROOT></ROOT>'))
        #self.cv = cascadingViewer(self.xmlWorking,self.schemaXML)
        #self.setCentralWidget(self.cv)

        #Cascading viewer
        self.cv = None
        #filename of xml that's currently in self.cv
        self.xmlFileName = None
        #xml that's in self.cv
        #self.xml = None
        #copy of what was initially in self.cv when the file was opened.
        #self.xmlWorking = None
        #Flag indicating whether the currently loaded xml has been changed
        #or not
        self.xmlChanged=None
        

    def initUI(self):
        """Initialise user interface."""

        self.setWindowTitle('Assessment Tool')

        self.statusBar().showMessage('Ready')
        
        #make a menu bar
        menuBar = self.menuBar()
        fileMenu = menuBar.addMenu('&File')
        dataMenu = menuBar.addMenu('&Data')
        reportMenu = menuBar.addMenu('&Plugins')

        
        #set up New action for the fileMenu
        #newAction = QtGui.QAction(QtGui.QIcon.fromTheme('filenew'),'New',self)
        newAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_FileIcon),'New',self)
        newAction.setShortcut('Ctrl+N')
        newAction.setStatusTip('New XML document')
        newAction.triggered.connect(self.newXML)
        
        #set up Open action for the fileMenu
        #openAction = QtGui.QAction(QtGui.QIcon.fromTheme('fileopen'),'Open',self)
        openAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DirIcon),'Open',self)
        openAction.setShortcut('Ctrl+O')
        openAction.setStatusTip('Open existing XML document')
        openAction.triggered.connect(self.openXML)
        
        #set up Save action for the fileMenu
        #saveAction = QtGui.QAction(QtGui.QIcon.fromTheme('filesave'),'Save',self)
        saveAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DialogSaveButton),'Save',self)
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save XML document')
        saveAction.triggered.connect(self.saveXML)
        
        #set up Save As action for the fileMenu
        saveAsAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DialogSaveButton),'Save As',self)
        saveAsAction.setShortcut('F12')
        saveAsAction.setStatusTip('Save XML document as ...')
        saveAsAction.triggered.connect(self.saveXMLAs)
        
        #set up Save As action for the fileMenu
        saveSelectionAsAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DialogApplyButton),'Save Visible',self)
        #saveSelectionAsAction.setShortcut('F12')
        saveSelectionAsAction.setStatusTip('Save XML for visible chemicals.')
        saveSelectionAsAction.triggered.connect(self.saveSelectionAs)

        #set up Close action for the fileMenu
        closeAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DialogCloseButton),'Close',self)
        closeAction.setShortcut('Ctrl+W')
        closeAction.setStatusTip('Close XML document')
        closeAction.triggered.connect(self.closeXML)

        aboutAction=QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_MessageBoxInformation),'About ...',self)
        aboutAction.triggered.connect(self.about)
        
        #set up exit action for the fileMenu
        exitAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_DialogCancelButton),'Exit',self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)
        
        #set up Download All action for the dataMenu
        downloadAllAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_ArrowDown), 'Download All',self)
        downloadAllAction.setStatusTip('Download all data in database.')
        if 'dbXml' in sys.modules:
            downloadAllAction.triggered.connect(self.downloadAll)

        #set up Upload xml action for the dataMenu
        uploadXmlAction = QtGui.QAction(self.style().standardIcon(\
                QtGui.QStyle.SP_ArrowUp), 'Upload xml',self)
        uploadXmlAction.setStatusTip('Upload new data in current xml to database.')
        if 'dbXml' in sys.modules:
            uploadXmlAction.triggered.connect(self.uploadXml)

        #set up Download List action for /the dataMenu
        #downloadListAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Download List',self)
        #downloadListAction.setStatusTip(\
        #        'Download data in database for a list of chemicals.')

        #set up Substructure action for the dataMenu
        #substructureAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Substructure search',self)
        #substructureAction.setStatusTip(\
        #        'Download data in database for chemicals having a particular substructure.')

        #set up Similarity action for the dataMenu
        #similarityAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Similarity search',self)
        #similarityAction.setStatusTip('Download data in database for '\
        #        'chemicals having a particular substructure.')

        ##set up Merge action for the dataMenu
        #mergeAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),'Merge',self)
        #mergeAction.setStatusTip('Merge XML from file into current list.')

        #set up Preassess action for the dataMenu
        #preassessAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Preassess',self)
        #preassessAction.setStatusTip('Run preassessment profiler on current '\
        #        'list of chemicals.')

        #Set up the reports. These are plugins (found in the plugin folder)
        print 'Initialising report plugins ...'


        self.reportActions = []
        for r in plugin.reportProvider.plugins:
            print r
            reportObject = r('sdfsdf')
            #print reportObject
            self.reportActions.append(QtGui.QAction(\
                    reportObject.returnReportType(),self))

            #reportMenu.addAction(self.reportActions[-1])

            del reportObject
        
            #make callback function
            callback = self.callbackFactory(r)
            #print callback()

            #need to make new variable r2 which is defined
            #within the lambda function so that r is stored.
            #This doesn't work for some reason.
            #callback = lambda r2=r: self.initiateReport(r2)


            #pass the name of the report to the slot
            self.reportActions[-1].triggered.connect(callback)

            #make callback object
            #cbObj = callbackObj(r,self.cv.model)
            #self.reportActions[-1].triggered.connect(cbObj.initiateReport)

        for reportAction in sorted(self.reportActions, key = lambda x: x.text()):
            reportMenu.addAction(reportAction)


        #r = plugin.reportProvider.plugins[0]
        #print '---'
        #for item in self.reportActions:
        #    print item.text()
        #print '-----'

        ##set up group report action for the reportMenu
        #groupReportAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Group report',self)
        #groupReportAction.setStatusTip('Generate grouping report for current '\
        #        'list of chemicals')
        #
        ##set up Tier I report action for the reportMenu
        #tierIReportAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Tier I report',self)
        #tierIReportAction.setStatusTip('Generate Tier I report for selected '\
        #        'chemical.')

        ##set up Tier I report action for the reportMenu
        #tierIReportAllAction = QtGui.QAction(QtGui.QIcon.fromTheme('down'),\
        #        'Tier I reports (All)',self)
        #tierIReportAllAction.setStatusTip('Generate Tier I report for all '\
        #        'chemicals.')

        fileMenu.addAction(newAction)
        fileMenu.addAction(openAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(saveAsAction)
        fileMenu.addAction(saveSelectionAsAction)
        fileMenu.addAction(closeAction)
        fileMenu.addAction(aboutAction)
        fileMenu.addAction(exitAction)

        dataMenu.addAction(downloadAllAction)
        dataMenu.addAction(uploadXmlAction)
        #dataMenu.addAction(downloadListAction)
        #dataMenu.addAction(substructureAction)
        #dataMenu.addAction(similarityAction)
        #dataMenu.addAction(mergeAction)
        #dataMenu.addAction(preassessAction)

        #add label to status bar
        self.visibleLabel = QtGui.QLabel()
        self.statusBar().addPermanentWidget(self.visibleLabel)

        #reportMenu.addAction(groupReportAction)
        #reportMenu.addAction(tierIReportAction)
        #reportMenu.addAction(tierIReportAllAction)

        self.setGeometry(50,50,800,700)
        #self.setGeometry(QtGui.QDesktopWidget().availableGeometry())

        #centre the window
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        self.show()

    def about(self):
        """Shows a messagebox containing commit number etc.
        
        Note that a special commit will have to be generated to change
        the commit number, since you don't know the commit number until
        you've committed it..."""
        QtGui.QMessageBox.information(self,'Cascading Viewer',\
            '''ECAT - the Electronic Chemical Assessment Tool.\n\n
Copyright The Department of the Environment,
Commonwealth of Australia, 2016.

Released under the GNU General Public License Version 2.''',\
            QtGui.QMessageBox.Ok)


    def closeXML(self):
        """Closes xml.
        
        Deletes the cascadingViewer widget and kills the autosave thread."""
        print 'close XML!'
        print self.cv

        #def closexml():
        #    #kill autosave thread
        #    print 'emitting stopAutoSave ...'
        #    self.cv.stopAutoSave.emit()

        #    self.cv.deleteLater()
        #    self.cv = None
        #    self.xmlFileName = None
        #    self.xmlChanged = None
        #    gc.collect()
        #    #objgraph.show_backrefs([self.xml],filename='xmlBackrefs.png')
        #    #del self.xml
        #    #self.xml = None

        #    return


        if self.cv:
            if self.xmlChanged == True:
                #prompt for save, then close.
                saveMsg = QtGui.QMessageBox()
                saveMsg.setText('The document '+self.xmlFileName+' has changed. Do you want to save your changes?')
                saveMsg.setStandardButtons(QtGui.QMessageBox.Save | \
                        QtGui.QMessageBox.Discard | \
                        QtGui.QMessageBox.Cancel)
                saveMsgOut = saveMsg.exec_()

                if saveMsgOut == QtGui.QMessageBox.Save:
                    self.saveXML()
                    #closexml()
                    self.cv.stopAutoSave.emit()
                    self.cleanUpCv()
                    return
                
                elif saveMsgOut == QtGui.QMessageBox.Discard:
                    #closexml()
                    self.cv.stopAutoSave.emit()
                    self.cleanUpCv()
                    return

                elif saveMsgOut == QtGui.QMessageBox.Cancel:
                    #do nothing
                    return 'cancelled'
                else:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                            'closeXML: Save went awry.',\
                        QtGui.QMessageBox.Ok)
                    return
            else:
                #closexml()
                self.cv.stopAutoSave.emit()
                self.cleanUpCv()
                return
        else:
            #Nothing to do, so do nothing
            return



    def saveXML(self):
        '''Saves XML.
        
        Always does the save, even if the xml has not changed. This is so
        That you can read in an unformatted (but valid) xml file, and 
        write a formatted one containing exactly the same data.'''
        def savexml(saveFileName):
            try:
                with open(saveFileName,'w') as saveFile:
                    sourceRoot = self.cv.model.\
                            createIndex(0,0,self.cv.model.rootItem).\
                            child(0,0).internalPointer().element
                    xml = copy.deepcopy(sourceRoot)

                    #print etree.tostring(tempXml,pretty_print=True)

                    print len(xml)
                    #print xml

                    saveFile.write('<ROOT>\n')
                    for item in xml:
                        a = etree.tostring(item,pretty_print=True)
                        saveFile.write(a)
                    saveFile.write('</ROOT>\n')

                    #Reset the xmlChanged flag.
                    self.xmlChanged = False
            except IOError:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                        'Couldn\'t save for some reason.',\
                        QtGui.QMessageBox.Ok)
                
               
        
        if not self.cv:
            return 
        else:
            if not self.xmlFileName:
                self.saveXMLAs()
                return 
            else:
                savexml(self.xmlFileName)
                return
            

    def saveXMLAs(self):
        """Save As functionality."""
        def savexml(saveFileName):
            try:
                with open(saveFileName,'w') as saveFile:
                    sourceRoot = self.cv.model.createIndex(0,0,self.cv.model.rootItem).\
                            child(0,0).internalPointer().element
                    xml = copy.deepcopy(sourceRoot)

                    #print etree.tostring(tempXml,pretty_print=True)

                    print len(xml)
                    #print xml

                    saveFile.write('<ROOT>\n')
                    for item in xml:
                        a = etree.tostring(item,pretty_print=True)
                        saveFile.write(a)
                    saveFile.write('</ROOT>\n')

                    #Reset the xmlChanged flag.
                    self.xmlChanged = False
                    self.xmlFileName = saveFileName
                    
            except IOError:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                        'Couldn\'t save for some reason.',\
                        QtGui.QMessageBox.Ok)
                

        if not self.cv:
            return
        else:
            saveFileName = QtGui.QFileDialog.getSaveFileName(self,'Save XML')
            if saveFileName != QtCore.QString(u''):
                savexml(saveFileName)
                print 'xml saved.'
                return
            else:
                #Need this for upload to know if cancel
                #was pressed.
                return QtCore.QString(u'')

            

    def saveSelectionAs(self):
        """Saves xml for the currently selected chemicals."""
        if self.cv:
            saveFileName = QtGui.QFileDialog.getSaveFileName(self,'Save XML')

            #bail if cancel button is clicked.
            if saveFileName == QtCore.QString(u''):
                return

            tempXml = etree.Element('ROOT')
            
            sourceRoot = self.cv.model.createIndex(0,0,self.cv.model.rootItem).\
                    child(0,0)
            proxyRoot = self.cv.chemProxyModel.mapFromSource(sourceRoot)
            print self.cv.model.rowCount(sourceRoot)
            print self.cv.chemProxyModel.rowCount(proxyRoot)

            for row in xrange(self.cv.chemProxyModel.rowCount(proxyRoot)):
                #print row
                sourceIdx = proxyRoot.model().mapToSource(proxyRoot.child(row,0))
                #print sourceIdx
                tempXmlChemical = copy.deepcopy(sourceIdx.internalPointer().element)
                tempXml.append(tempXmlChemical)

            #print etree.tostring(tempXml,pretty_print=True)
            


            try:
                with open(saveFileName,'w') as saveFile:
                    tempXmlString = etree.tostring(tempXml,pretty_print=True)
                    saveFile.write(tempXmlString)


            except IOError:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                        'Couldn\'t save for some reason.',\
                        QtGui.QMessageBox.Ok)
                    return

            reply = QtGui.QMessageBox.question(self,'Message',
                    'Selection saved as xml. Would you like to open the file you just saved?',QtGui.QMessageBox.Yes|QtGui.QMessageBox.No,QtGui.QMessageBox.Yes)
            if reply == QtGui.QMessageBox.Yes:
                if self.xmlChanged == True:
                    reply2 = QtGui.QMessageBox.question(self,'Message',
                            self.xmlFileName+' has changed. Would you like to save your changes before closing it?',QtGui.QMessageBox.Yes| QtGui.QMessageBox.No,QtGui.QMessageBox.No)
                    if reply2 == QtGui.QMessageBox.Yes:
                        self.saveXMLAs()
                    elif reply2 ==QtGui.QMessageBox.No:
                        pass
                    else:
                        print 'Should be QtGui.QMessageBox.Yes or QtGui.QMessageBox.No.'
                        return

                #close and stop autosave.
                print 'emitting stopAutoSave ...'

                self.cv.stopAutoSave.emit()
                self.cleanUpCv()
                #self.cv.deleteLater()
                #self.cv = None
                #self.xmlFileName = None
                #self.xmlChanged = None
                ##self.xml = None
                #gc.collect()

                self.openXML(saveFileName)
            else:
                return


        else:
            return

    def openXML(self,fname=None):
        '''Opens XML database.
        
        This borrows the QtGui.QMessageBox standard button enums.
        Creates a cascadingViewer widget (which kicks off the autosave thread). 
        '''

        #pprint(objgraph.typestats())

        def getFnameIfNecessary(fname):
            if fname == None or fname == False:
                #get filename to open. Otherwise, we already have it.
                fname = QtGui.QFileDialog.getOpenFileName(self,'Open File','/home/dbrittain/projects/pyqt/cascading_viewer')
            return fname
        
        def myOpenXml(fname):

            try:
                parser = etree.XMLParser(remove_blank_text=True)
                tempXmlTree = etree.parse(str(fname),parser)

            except:
                QtGui.QMessageBox.information(self,'Cascading Viewer',\
                    'The file you tried to open does not exist \
                     or is not an XML document.',\
                    QtGui.QMessageBox.Ok)
                return 'cancelled'
            
            #At this point, we have valid xml. Check against schema.
            tempXml = tempXmlTree.getroot()

            #objgraph.show_growth()
            #testing

            if self.schemaRNG.validate(tempXml):
                print 'xml validates'
                #self.xml = tempXml
                #self.xmlWorking = deepcopy(self.xml)
                self.xmlFileName = fname


                self.cv = cascadingViewer(tempXml,self.schemaXML,self.xmlFileName,self)
                self.xmlChanged = False
                self.setCentralWidget(self.cv)

                #connect dataChanged signal so that xmlChanged
                #is set to true as soon as data is altered.
                self.cv.model.dataChanged.connect(self.setXmlChangedTrue)
                
                self.cv.filterApplied.connect(self.changeNumberVisible)
                self.changeNumberVisible()
                #objgraph.show_backrefs([tempXml],filename='xmlBackrefs.png')

                return
            else:
                QtGui.QMessageBox.information(self,'Cascading Viewer',\
                    'The XML file you tried to open does not validate against '\
                    'the schema.',\
                    QtGui.QMessageBox.Ok)
                return

        #OpenXML proper
        if self.cv:
            if self.xmlChanged == True:
                #prompt for save, then close.
                saveMsg = QtGui.QMessageBox()
                saveMsg.setText('The document '+self.xmlFileName+' has changed. Do you want to save your changes?')
                saveMsg.setStandardButtons(QtGui.QMessageBox.Save | \
                        QtGui.QMessageBox.Discard | \
                        QtGui.QMessageBox.Cancel)
                saveMsgOut = saveMsg.exec_()

                if saveMsgOut == QtGui.QMessageBox.Save:
                    self.saveXML()
                    self.cv.stopAutoSave.emit()
                    self.cleanUpCv()
                elif saveMsgOut == QtGui.QMessageBox.Discard:
                    self.cv.stopAutoSave.emit()
                    self.cleanUpCv()
                elif saveMsgOut == QtGui.QMessageBox.Cancel:
                    #do nothing
                    return 'cancelled'
                else:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                            'Save went awry.',\
                        QtGui.QMessageBox.Ok)
                    return
        fname = getFnameIfNecessary(fname)
        if fname != QtCore.QString(u''):
            myOpenXml(fname)
        return

    def newXML(self):
        print 'New xml!'
        if self.cv:
            self.closeXML()

        with open(os.path.join(os.getcwd(),'autosave/new.xml'),'w') as f:
            f.write('<ROOT>\n</ROOT>')

        self.openXML(os.path.join(os.getcwd(),'autosave/new.xml'))


    #Reimplement closeEvent() event handler.
    def closeEvent(self,event):
        """Re-implemented closeEvent handler."""
        print 'close XML on close event!'
        #The event blocks signals. Just get rid of it.
        event.ignore()

        if self.cv:
            if self.xmlChanged == True:
                #prompt for save, then close.
                saveMsg = QtGui.QMessageBox()
                saveMsg.setText('The document '+self.xmlFileName+' has changed. Do you want to save your changes?')
                saveMsg.setStandardButtons(QtGui.QMessageBox.Save | \
                        QtGui.QMessageBox.Discard | \
                        QtGui.QMessageBox.Cancel)
                saveMsgOut = saveMsg.exec_()

                if saveMsgOut == QtGui.QMessageBox.Save:
                    reply =self.saveXMLAs()
                    if reply == QtCore.QString(u''):
                        return
                    else:
                        self.cv.destroyed.connect(QtGui.QApplication.quit)
                        self.cv.stopAutoSave.emit()
                        self.cleanUpCv()
                        #self.cv.autoSaveThread.quit()
                        #self.cv.autoSaveThread.wait()
                        #self.cv.autoSaveThread.deleteLater()
                        #QtGui.QApplication.quit()
                
                elif saveMsgOut == QtGui.QMessageBox.Discard:
                    #closexml()
                    self.cv.destroyed.connect(QtGui.QApplication.quit)
                    self.cv.stopAutoSave.emit()
                    self.cleanUpCv()
                    #self.cv.autoSaveThread.quit()
                    #self.cv.autoSaveThread.wait()
                    #self.cv.autoSaveThread.deleteLater()
                    #QtGui.QApplication.quit()
                    return

                elif saveMsgOut == QtGui.QMessageBox.Cancel:
                    #do nothing
                    return
                else:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                            'closeXML: Save went awry.',\
                        QtGui.QMessageBox.Ok)
                    return
            else:
                #closexml()
                self.cv.destroyed.connect(QtGui.QApplication.quit)
                self.cv.stopAutoSave.emit()
                self.cleanUpCv()
                #self.cv.autoSaveThread.quit()
                #self.cv.autoSaveThread.wait()
                #self.cv.autoSaveThread.deleteLater()
                #QtGui.QApplication.quit()

                

        else:
            #Nothing to do, so do nothing, but accept event
            QtGui.QApplication.quit()

    def callbackFactory(self,r):
        """Callback factory for plugins."""
        return lambda: self.initiateReport(r)

    def initiateReport(self,reportClass):
        """Runs a plugin."""
        if self.cv:
            print reportClass
            selectedChem = self.cv.chemSelectionModel.currentIndex()
            selectedChem = self.cv.chemProxyModel.mapToSource(selectedChem)

            selectedDp = self.cv.dataPointSelectionModel.currentIndex()
            selectedDp = self.cv.dataPointProxyModel.mapToSource(selectedDp)
            
            selectedQ = self.cv.qualifierSelectionModel.currentIndex()
            selectedQ = self.cv.qualifierProxyModel.mapToSource(selectedQ)

            print self.cv.model,selectedChem,selectedDp,selectedQ,self
            report = reportClass((self.cv.model,selectedChem,selectedDp,selectedQ,self))
            #report.makeReport()
            report.displayReport()

    def uploadXml(self):
        """Uploads xml to database. 
        
        The class that does all the work is found in dbXml.py."""
        #This is ugly, but it works.
        if self.cv != None:

            #Do some basic checks before we do anything else.
            print 'checking VARHCAR2 lengths ...'
            
            xpathStrings = ['//DP_QUAL_VALUE1/text()','//DP_QUAL_VALUE2/text()',\
                            '//Q_QUAL_VALUE1/text()','//Q_QUAL_VALUE2/text()']
            for xpathString in xpathStrings:
                qualValStrings = self.cv.model.rootItem.element[0].xpath(xpathString)
                for item in qualValStrings:
                    if len(item) > 4000:
                        printString = '\"'+str(item) + '\" is longer than 4000 characters.'
                        reply = QtGui.QMessageBox.critical(self,'Danger, Will Robinson!',
                            printString,QtGui.QMessageBox.Ok,QtGui.QMessageBox.Ok)
                        return


            reply = QtGui.QMessageBox.question(self,'Danger, Will Robinson!',
                    'Are you really sure you want to upload the contents of this xml document to the database?',QtGui.QMessageBox.Yes|
                    QtGui.QMessageBox.No,QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.No:
                return


            saveReturn = self.saveXMLAs()
            if saveReturn == QtCore.QString(u''):
                #Bail if you hit cancel.
                return
            #close, but save filename so we can read it.
            fname = str(self.xmlFileName)
            
            #self.cv.deleteLater()
            print 'emitting stopAutoSave ...'
            self.cv.stopAutoSave.emit()
            self.cleanUpCv()

            #self.cv.deleteLater()
            #self.cv = None
            #self.xmlFileName = None
            #self.xmlChanged = None
            ##self.xml = None
            #gc.collect()




            #get output name for zipped xml
            outFName = fname
            i = 1
            while os.path.exists(outFName+'.bz2'):
                outFName = fname +'.'+str(i)
                i=i+1
            
            print 'Compressing data ...'
            with open(fname,'rb') as ip:
                with bz2.BZ2File(outFName+'.bz2','wb',compresslevel=9) as op:
                    copyfileobj(ip,op)

            print 'Backed up xml to '+outFName+'.bz2'

            
            
            #finally, do the upload.
            #check to see that the file exists
            if os.path.exists(outFName+'.bz2'):
                xml = etree.parse(fname)
                xml = xml.getroot()
                up = dbXml.uploadfromxml(self.schemaRNG,self)
                up.upload(xml)
                os.remove(fname)
                #import cx_Oracle
                #con = cx_Oracle.connect('cast_owner','cast_owner','intltest')
                #exit()

            


    def downloadAll(self):
        """Sucks down all of the data from the database and save in xml format."""
        filename = QtGui.QFileDialog.getSaveFileName(self,'Download xml to:')
        if filename == QtCore.QString(u''):
            return
        with dbXml.downloadToXml() as dw:
            #dw = dbXml.downloadToXml()
            dw.pullXml(filename)



        reply = QtGui.QMessageBox.question(self,'Message',
                'Xml downloaded. Would you like to open it?',QtGui.QMessageBox.Yes|
                QtGui.QMessageBox.No,QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            self.openXML(filename)
            

    def setXmlChangedTrue(self):
        """sets self.xmlChanged to True."""
        print 'setXmlChangedTrue'
        self.xmlChanged = True

    def changeNumberVisible(self):
        """Changes the count of visible chemicals in the status bar."""
        print 'changeNumberVisible'
        cpm = self.cv.chemProxyModel
        idx = cpm.index(0,0,QtCore.QModelIndex())
        print idx == QtCore.QModelIndex()
        print 'cpm.rowCount(idx): ',cpm.rowCount(idx)
        sourceIdx = cpm.mapToSource(idx)
        if cpm.rowCount(idx) == 0:
            message = 'No chemicals to display.'
        elif cpm.rowCount(idx) == 1:
            message = str(cpm.rowCount(idx))+ ' chemical of '+str(sourceIdx.model().rowCount(sourceIdx))+' chemicals visible'
        else:
            message = str(cpm.rowCount(idx)) +' chemicals of '+str(sourceIdx.model().rowCount(sourceIdx))+' chemicals visible'
            

        #message = 'Hello!'
        #print message
        #self.statusBar().showMessage(message)
        self.visibleLabel.setText(message)

        #also, if there is no chemical selected, select the first one.
        #print 'Has selection?', self.cv.chemSelectionModel.hasSelection()
        if self.cv.chemSelectionModel.hasSelection() == False:
            self.cv.chemSelectionModel.setCurrentIndex(self.cv.chemProxyModel.mapFromSource(self.cv.model.createIndex(0,0,self.cv.model.rootItem)).child(0,0),QtGui.QItemSelectionModel.SelectCurrent)
    
    def cleanUpCv(self):
        """This slot is connected to a cv's destroyed signal."""
        print 'Cleaning up cv ...'
        self.cv.deleteLater()
        self.cv = None
        self.xmlFileName = None
        self.xmlChanged = None
        gc.collect()
        return

class newDoubleValidator(QtGui.QDoubleValidator):
    """Modified QDoubleValidator.
    
    There is a bug in qt such that QDoubleValidator doesn't take any notice
    of RejectGroupSeparator in the current locale. This means it always accepts
    strings of the form "2.3,4.5,3333.444" because it thinks it's a list of doubles.

    This has been fixed in Qt version 5.4.1, see https://bugreports.qt.io/browse/QTBUG-42522 for details.

    More precisely:
    newLocale = QtCore.QLocale(QtCore.QLocale.English,QtCore.QLocale.Australia)
    newLocale.setNumberOptions(newLocale.numberOptions()|QtCore.QLocale.RejectGroupSeparator)
    QtCore.QLocale.setDefault(newLocale)
    
    has no effect on QDoubleValidator.
    
    This class just checks to see if there's a comma in the input string before
    running the parent class' validate function."""
    def __init__(self):
        super(newDoubleValidator,self).__init__()

    def validate(self,qstr,pos):
        """Re-implemented validate."""
        if ',' in str(qstr):
            return (QtGui.QValidator.Invalid,pos)
        else:
            return super(newDoubleValidator,self).validate(qstr,pos)


def main():
    """Kicks off the event loop."""
    app = QtGui.QApplication(sys.argv)
    #cv = cascadingViewer()


    at = assessmentTool()
    
    sys.exit(app.exec_())
    
if __name__ == '__main__':

    main()
