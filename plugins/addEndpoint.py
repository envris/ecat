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

class endpointWidget(QtGui.QWidget):
    """Endpoint widget."""
    def __init__(self):
        #call constructor of parent class
        super(endpointWidget,self).__init__()
        
        self.columnHeaderDict = self.makeColumnHeaderDict()

        #call method that generates the UI.
        self.initUI()
    
    def initUI(self):
        """Initialises user interface."""


        self.valueLineEdit = QtGui.QLineEdit()
        dv = QtGui.QDoubleValidator()
        self.valueLineEdit.setValidator(dv)

        self.unit = QtGui.QComboBox()
        self.unit.addItem('mg/L')
        self.unit.addItem('mg/kg')
        self.unit.addItem('ug/kg bodyweight')
        self.unit.addItem('ug/kg soil')
        self.unit.addItem('ug/bee')
        self.unit.addItem('ug/organism')
        self.unit.addItem('kg/hectare')

        self.trophicCombo = QtGui.QComboBox()
        self.trophicCombo.setEditable(True)
        self.trophicCombo.addItems(['Fish','Invertebrate','Algae'])
        self.trophicCombo.currentIndexChanged[str].connect(self.onTrophicIndexChanged)

        self.endpointNameCombo = QtGui.QComboBox()
        self.endpointNameCombo.addItems(['LC50','EC50'])
        self.endpointNameCombo.setEditable(True)

        self.durationLineEdit = QtGui.QLineEdit()
        self.durationLineEdit.setValidator(dv)

        self.durationUnitCombo = QtGui.QComboBox()
        self.durationUnitCombo.addItem('hour')
        self.durationUnitCombo.addItem('day')
        self.durationUnitCombo.addItem('week')
        self.durationUnitCombo.addItem('month')
        self.durationUnitCombo.addItem('year')

        self.pivotalCombo = QtGui.QComboBox()
        self.pivotalCombo.addItems(['N','Y'])

        self.trophicDict = {'Fish':['Pimephales promelas','Oncorhynchus mykiss',\
                    'Oryzias latipes'],\
                    'Invertebrate':['Daphnia magna'],\
                    'Algae':['Pseudokirchneriella subcapitata',\
                    'Chlorella pyrenoidosa']}
        self.speciesList = []
        for tl in self.trophicDict.keys():
            for item in self.trophicDict[tl]:
                self.speciesList.append(item)
        self.speciesCombo = QtGui.QComboBox()

        self.speciesCombo.setEditable(True)

        self.acuteChronicCombo = QtGui.QComboBox()
        self.acuteChronicCombo.addItem('Acute')
        self.acuteChronicCombo.addItem('Chronic')

        self.measuredCombo = QtGui.QComboBox()
        self.measuredCombo.addItem('Measured')
        self.measuredCombo.addItem('Calculated')
        self.measuredCombo.addItem('Read-across')
        self.measuredCombo.addItem('Analogue')

        #initalize species combo
        self.onTrophicIndexChanged(self.trophicCombo.currentText())




        #formLayout = QtGui.QFormLayout()
        #formLayout.addRow('Endpoint value (mg/L)',self.valueLineEdit)
        #formLayout.addRow('Trophic Level',self.trophicCombo)
        #formLayout.addRow('Endpoint',self.endpointNameCombo)
        #formLayout.addRow('SAR',self.SARLineEdit)
        #formLayout.addRow('Duration (hours)',self.durationLineEdit)
        #formLayout.addRow('Pivotal endpoint?',self.pivotalCombo)

        gl = QtGui.QGridLayout()
        gl.addWidget(QtGui.QLabel('Value:'),0,0)
        gl.addWidget(self.valueLineEdit,0,1)
        gl.addWidget(self.unit,0,2)

        gl.addWidget(self.durationLineEdit,1,0)
        gl.addWidget(self.durationUnitCombo,1,1)
        gl.addWidget(self.endpointNameCombo,1,2)
        
        gl.addWidget(QtGui.QLabel('Organism type:'),2,0)
        gl.addWidget(self.trophicCombo,2,1)

        gl.addWidget(QtGui.QLabel('Species:'),3,0)
        gl.addWidget(self.speciesCombo,3,1)

        gl.addWidget(QtGui.QLabel('Acute or Chronic?'),4,0)
        gl.addWidget(self.acuteChronicCombo,4,1)

        gl.addWidget(QtGui.QLabel('Measured?'),5,0)
        gl.addWidget(self.measuredCombo,5,1)

        gl.addWidget(QtGui.QLabel('Pivotal Endpoint?'),6,0)
        gl.addWidget(self.pivotalCombo,6,1)



        #initialize data


        #self.setLayout(formLayout)
        self.setLayout(gl)

        self.setGeometry(300,300,350,300)
        #self.move(300,150)
        self.setWindowTitle('Review')
        self.show()

    def onTrophicIndexChanged(self,trophicString):
        "Alter species combo depending on trophic level selection."

        initialText =  str(self.speciesCombo.currentText())
        ts = str(trophicString)


        self.speciesCombo.clear()
        if initialText in self.speciesList or initialText == '':
            if ts in self.trophicDict.keys():
                #set contents of species combo
                for item in self.trophicDict[ts]:
                    self.speciesCombo.addItem(item)
        else:
            self.speciesCombo.addItem(initialText)






    def makeQYN(self,name,YorN):
        """Make Y/N qualifier."""
        if YorN not in set(['Y','N']):
           print 'makeQYN: YorN must be \'Y\' or \'N\''
           exit()
        q = etree.Element(name)

        qid =  etree.SubElement(q,'QUALIFIER_ID')
        qid.text = 'AUTO'

        dpid =  etree.SubElement(q,'DATA_POINT_ID')
        dpid.text = 'AUTO'

        kc = etree.SubElement(q,'Q_L1_KINGDOM_CODE')
        kc.text = 'Q_YN'

        cb = etree.SubElement(q,'CREATED_BY')
        cb.text = 'AUTO'

        cd = etree.SubElement(q,'CREATED_DATE')
        cd.text = 'AUTO'

        confidential = etree.SubElement(q,'Q_CONFIDENTIAL1_YN')
        confidential.text = 'N'

        qv1 = etree.SubElement(q,'Q_QUAL_VALUE1')
        qv1.text = YorN

        return q

    def makeQVarchar2(self,dpTag,dpData):
        """Make VARCHAR2 qualifier."""
        #dpType = self.columnHeaderDict[colHeader]['type']
        #dpTag = 'SAR'
        #print colHeader, dpData,dpTag
        errorReturn = None

        if dpData != '' and dpData != '#N/A' and dpData != None:

            Xml = etree.Element(dpTag)

            dpid =  etree.SubElement(Xml,'QUALIFIER_ID')
            dpid.text = 'AUTO'

            chemid = etree.SubElement(Xml,'DATA_POINT_ID')
            chemid.text = 'AUTO'

            l1 = etree.SubElement(Xml,'Q_L1_KINGDOM_CODE')
            l1.text = 'Q_VARCHAR2'

            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'

            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'

            confidential = etree.SubElement(Xml,'Q_CONFIDENTIAL1_YN')
            confidential.text = 'N'

            varchar2 = etree.SubElement(Xml,'Q_QUAL_VALUE1')
            varchar2.text = dpData

            return Xml
        else:
            return None

    def makeQFloat(self,qTag,unit,data):
        """Make float qualifier."""
        #Only makes single floats (used for 96-hr in endpoints)
        if data != '' and data != '#N/A' and data != None:
            Xml = etree.Element(qTag)

            dpid =  etree.SubElement(Xml,'QUALIFIER_ID')
            dpid.text = 'AUTO'

            chemid = etree.SubElement(Xml,'DATA_POINT_ID')
            chemid.text = 'AUTO'

            l1 = etree.SubElement(Xml,'Q_L1_KINGDOM_CODE')
            l1.text = 'Q_FLOAT'

            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'

            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'

            qqv1 = etree.SubElement(Xml,'Q_QUANT_VALUE1')
            qqv1.text = data

            confidential = etree.SubElement(Xml,'Q_CONFIDENTIAL1_YN')
            confidential.text = 'N'

            qu1 = etree.SubElement(Xml,'Q_UNIT_CODE1')
            qu1.text = unit

            return Xml
        else:
            return None



    def makeFloat(self,colHeader,unit,dpData):
        """Make float data point."""


        dpType = self.columnHeaderDict[colHeader]['type']
        dpTag = self.columnHeaderDict[colHeader]['tag']
        #print colHeader, dpData,dpTag
        errorReturn = None

        if dpData != '' and dpData != '#N/A' and dpData != None:

            Xml = etree.Element(dpTag)
            
            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'
            
            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'
           

            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_FLOAT'
                
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'
    
            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
            
            #special cases:
            if colHeader == 'AHVICL Threshold Range (2006)':
                if dpData == '1000 to 9999 tonnes':
            
                    f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                    f1.text = '1000.0'
                    
                    f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                    f2.text = '9999.0'

                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                    c1.text = 'Y'
                    
                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                    c1.text = 'Y'

                    
                    q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                    q1.text = '>='
            
                    q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                    q2.text = '<='

                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = 'tonne'
            
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = 'tonne'
                    
            
                elif dpData == '10,000 to 99,999 tonnes':
            
                    f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                    f1.text = '10000.0'
                    
                    f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                    f2.text = '99999.0'
                    

                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                    c1.text = 'Y'
                    
                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                    c1.text = 'Y'
                    
                    q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                    q1.text = '>='

                    q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                    q2.text = '<='
            
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = 'tonne'
            
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = 'tonne'
                    
            
                elif dpData == '100,000 to 1,000,000 tonnes':
                    f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                    f1.text = '100000.0'
                    
                    f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                    f2.text = '999999.0'
                    
                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                    c1.text = 'Y'
                    
                    c1 =  etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                    c1.text = 'Y'

                    q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                    q1.text = '>='
            
                    q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                    q2.text = '<='


                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = 'tonne'
            
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = 'tonne'
                    
            
                elif dpData == '' or dpData == 'Not on list':
                    return None
                else:
                    #errorReturn = [self.data[index][self.fielddict['CAS-RN']],\
                    #        colHeader,dpData]
                    errorReturn = [colHeader,dpData]


            
            
            #parse
            elif re.match('^ *(?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?$',dpData):
                #it's a straight number.

                f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                f1.text = dpData


                confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                confidential.text = 'N'
    
                if unit != None:
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = unit
                #if re.match('^(?:<|<=|>|>=)? *(?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?$',dp):
                #pass
            elif re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData):
                #(>,>=,<,<=)
                gtlt = re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData).group(1)
                num = re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData).group(2)
            
                #lower bound.
                if gtlt == '>' or gtlt == '>=':
                    f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                    f1.text = num
                    
    
                    confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                    confidential.text = 'N'
        
                    q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                    q1.text = gtlt
            
                    #units
                    if unit != None:
                        u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                        u1.text = unit
                

                #upper bound
                elif gtlt == '<' or gtlt == '<=':
                    f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                    f2.text = num
                    
                    confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                    confidential.text = 'N'
        
                    
                    q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                    q2.text = gtlt
                    
                    #units
                    if unit != None:
                        u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                        u2.text = unit
            
            elif re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?) *to *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData):
                # a to b

                atobmatch = re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?) *to *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData)                


                f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                f1.text = atobmatch.group(1)

                f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                f2.text = atobmatch.group(2)


                confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                confidential.text = 'N'
                
                confidential2 = etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                confidential2.text = 'N'
    
            
                q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                q1.text = '>='

                q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                q2.text = '<='
                
                if unit != None:
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = unit
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = unit

            elif re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?) *- *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData):
                # a - b

                atobmatch = re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?) *- *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)?[0-9]+)?)$',dpData)


                f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                f1.text = atobmatch.group(1)

                f2 = etree.SubElement(Xml,'DP_QUANT_VALUE2')
                f2.text = atobmatch.group(2)



                confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                confidential.text = 'N'
                
                confidential2 = etree.SubElement(Xml,'DP_CONFIDENTIAL2_YN')
                confidential2.text = 'N'
    
                q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                q1.text = '>='
                
                q2 = etree.SubElement(Xml,'DP_QUAL_VALUE2')
                q2.text = '<='

                if unit != None:
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = unit
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = unit

            
            elif dpData == '':
                return None

            else:
                #errorReturn = [self.data[index][self.fielddict['CAS-RN']],\
                #    colHeader,dpData]
                errorReturn = [colHeader,dpData]

            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)


            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None



    def makeFloatException(self,colHeader,dpData):
        """Make exception-type float data point (e.g. "Miscible" for water solubility)."""


        dpType = self.columnHeaderDict[colHeader]['type']
        dpTag = self.columnHeaderDict[colHeader]['tag']
        #print colHeader, dpData,dpTag
        errorReturn = None

        if dpData != '' and dpData != '#N/A' and dpData != None:

            Xml = etree.Element(dpTag)
            
            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'
            
            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'
           

            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_FLOAT'
                
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'
    
            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
            
            
            if dpData in ['No effects at saturation','NES',\
                    'No effect at saturation', 'No effect at saturation.',\
                    'No Effects at Saturation']:
                confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                confidential.text = 'N'

                q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                q1.text = 'No effects at saturation'


            else:
                #errorReturn = [self.data[index][self.fielddict['CAS-RN']],\
                #    colHeader,dpData]
                errorReturn = [colHeader,dpData]

            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)


            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None

    
    def makeColumnHeaderDict(self):
        """Makes dictionary of column headers."""

        wbInfo = [{'CAS-RN': {'type': 'DP_INTEGER', 'unit': None, 'tag' : 'CASNUMBER'},
        'RQ': {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'RQ'},
        'Release Mitigation Factor' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'RMF'},
        'Volume (T) (8 May 2013)' : {'type' : 'DP_FLOAT', 'unit' : 'tonne', 'tag' : 'Volume'},
        'Mitigated PEC(river) (ug/L)' : {'type' : 'DP_FLOAT', 'unit' : 'ug/L', 'tag' : 'MitigatedPEC'},
        'PNEC (ug/L)': {'type' : 'DP_FLOAT', 'unit' : 'ug/L', 'tag' : 'PNEC' },
        'Persistent?' : {'type' : 'DP_VARCHAR2CHOICE', 'unit' : None, 'tag' : 'P','choices' : ['Y','N','P*','Uncertain'] },
        'Bioaccumulative?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'B' },
        'Toxic?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'T' }},
        
        {'Assessment conclusion' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'Conclusion' },
        'Comments' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'Comments' },
        'Internal Notes (not for publication)' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'InternalNotes' },
        'Group' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'Group' },
        'Assessment Status' : {'type' : 'DP_VARCHAR2CHOICE', 'unit' : None, 'tag' : 'Status','choices' : ['Assessment - T1', 'Complete - T1','Excluded - T1', 'Grouped - T1','Grouped - T2' ,'Initial Review - T1', 'Peer Review - T1', 'Peer Review - T2', 'Published - T1', 'Published - T2', 'Ungrouped - T1', 'Ungrouped - T2', 'Ungrouped - T3'] },
        'Batch no.' : {'type' : 'DP_INTEGER', 'unit' : None, 'tag' : 'BatchNo' },
        'Assessor' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'Assessor' }},
        
        {'Reasonable worst case exposure scenario appropriate?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'ReasonableWorstCaseExposureScenarioAppropriate' },
        #'Number of HVICL Uses' : {'type' : 'DP_INTEGER', 'unit' : None, 'tag' : 'numberOfHVICLUses'},
        'AHVICL Uses (2006)' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'HVICLUse' },
        'AHVICL Threshold Range (2006)' : {'type' : 'DP_FLOAT', 'unit' : 'tonne', 'tag' : 'HVICLRange' },
        #'SPIN UC62 - Total Tonnage' : {'type' : 'DP_FLOAT', 'unit' : 'tonne', 'tag' : 'SPINUC62' },
        #'Top Three SPIN Uses (UC62)' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'SPINUsesUC62'},
        #'Number of SPIN Uses (UC62)' : {'type' : 'DP_INTEGER', 'unit': None, 'tag' : 'NumberOfSpinUsesUC62'},
        #'Number of detailed SPIN Uses' : {'type' : 'DP_INTEGER', 'unit': None, 'tag' : 'NumberOfDetailedSpinUses'},
        'SPIN Uses (UC62)' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'SPINUsesUC62'},
        'Detailed SPIN Uses' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'DetailedSPINUses'}} ,
        {'NICNAS ID' : {'type' : 'DP_INTEGER', 'unit' : None, 'tag' : 'NICNAS_ID' },
        'AICS Name': {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'AICSName' },
        'Common name' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'CommonName' },
        'IMAP Status' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'IMAPStatus' },
        'Chemical class' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'ChemicalClass' },
        'SMILES String' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'SMILES' },
        'Molecular formula' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'MolFormula' } ,
        'Molecular weight (g/mol)' : {'type' : 'DP_FLOAT', 'unit' : 'g/mol', 'tag' : 'MolWeight' }},
        {'High concern use for environment?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'highConcernUseForEnvironment' },
        'Perfluorinated? (NICNAS master list 12 Oct 2012)' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'Perfluorinated' },
        'Montreal' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'Montreal' },
        'SGG' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'SGG' },
        'Rotterdam' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'Rotterdam' },
        'Stockholm' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'Stockholm' },
        'REACH (SVHCs)' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'AnnexeXIVREACH' },
        'EDC (US EPA)' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'EDCUSEPA' },
        'EDC (Europe)' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'EDCEurope' },
        'On DSL' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'OnDSL' },
        'DSL P' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'DSL_P' },
        'DSL B' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'DSL_B' },
        'DSL iT' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'DSL_iT' }},
        {'log Kow' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'logKow' },
            'Water solubility  (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'WaterSolubility','exceptions':['Miscible'] },
        'Melting point (deg C)' : {'type' : 'DP_FLOAT', 'unit' : 'deg C', 'tag' : 'MeltingPoint' },
        'Vapour pressure (Pa)' : {'type' : 'DP_FLOAT', 'unit' : 'Pa', 'tag' : 'VapourPressure' },
        'pKa' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'pKa' },
        'Ionisable in the environment?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'IonisableInTheEnvironment'},
        'Phys Chem Notes' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'PhysChemNotes' }},
        {'Persistent?' : {'type' : 'DP_VARCHAR2CHOICE', 'unit' : None, 'tag' : 'P','choices' : ['Y','N','P*','Uncertain']},
        'Obs. BOD (301C)' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'ObsBOD301C' },
        'BOD (301C)' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'CalcBOD301CCatalogic' },
        'Primary halflife in water (days; Catalogic 301C)' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'PrimaryHalfLifeWater' },
        'Reasons for Categorisation (P)' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'ReasonsForPCategorisation' },
        'Bioaccumulative?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'B' },
        'Reasons for Categorisation (B)' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'ReasonsForBCategorisation' },
        'Toxic?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'T' },
        'Toxicity Notes' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'ToxicityNotes' }},
        {'Fish ECOSAR (mg/L; Neutral Organics SAR)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARFishNeutralOrganics','exceptions': ['No effects at saturation'] },
                'Daphnia ECOSAR (mg/L; Neutral Organics SAR)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARDaphniaNeutralOrganics','exceptions': ['No effects at saturation'] },
                'Algae ECOSAR (mg/L; Neutral Organics SAR)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARAlgaeNeutralOrganics','exceptions': ['No effects at saturation'] } },
        {'ECOSAR Acute fish endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARAcuteFishEndpoint','exceptions': ['No effects at saturation'] },
        'Fish SAR' : {'type' : 'Q_VARCHAR2', 'unit' : None, 'tag' : 'ECOSARAcuteFishSAR' },
        'ECOSAR Acute invertebrate endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARAcuteInvertebrateEndpoint','exceptions': ['No effects at saturation'] },
        'Invertebrate SAR' : {'type' : 'Q_VARCHAR2', 'unit' : None, 'tag' : 'ECOSARAcuteInvertebrateSAR' },
        'ECOSAR Acute algae endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'ECOSARAcuteAlgaeEndpoint','exceptions': ['No effects at saturation'] },
        'Algae SAR' : {'type' : 'Q_VARCHAR2', 'unit' : None, 'tag' : 'ECOSARAcuteAlgaeSAR' }},
        {'Fish endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'measuredAcuteFishEndpoint','exceptions': ['No effects at saturation'] },
                'Invertebrate endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'measuredAcuteInvertebrateEndpoint','exceptions': ['No effects at saturation'] },
                'Algae endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'measuredAcuteAlgaeEndpoint','exceptions': ['No effects at saturation'] },
        'Toxicity Notes' : {'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'ToxicityNotes' }},
        
        {'Pivotal endpoint (mg/L)' : {'type' : 'DP_FLOAT', 'unit' : 'mg/L', 'tag' : 'PivotalEndpoint','exceptions': ['No effects at saturation'] },
        'Pivotal endpoint type' : {'type' : 'DP_VARCHAR2', 'unit' : None, 'tag' : 'PivotalEndpointType' },
        'Assessment factor' :{'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'AssessmentFactor' } ,
        'AF Notes' :{'type' : 'DP_CLOB', 'unit' : None, 'tag' : 'AssessmentFactorNotes' } ,
        'PNEC (ug/L)' :{'type' : 'DP_FLOAT', 'unit' : 'ug/L', 'tag' : 'PNEC' } },
        
        {'Default volume? (8 May 2013)' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'DefaultVolume' },
        'Volume (T) (8 May 2013)' :{'type' : 'DP_FLOAT', 'unit' : 'tonne', 'tag' : 'Volume' },
        'PEC(river) (ug/L ; 8 May 2013)' : {'type' : 'DP_FLOAT', 'unit' : 'ug/L', 'tag' : 'PECRiver' },
        'Release Mitigation Factor' :{'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'RMF' } ,
        'Mitigated PEC(river) (ug/L)' :{'type' : 'DP_FLOAT', 'unit' : 'ug/L', 'tag' : 'MitigatedPEC' } },
        {'RQ' :{'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'RQ' } ,
        'Release Mitigation Factor' : {'type' : 'DP_FLOAT', 'unit' : None, 'tag' : 'RMF' },
        'Volume (T) (8 May 2013)' : {'type' : 'DP_FLOAT', 'unit' : 'tonne', 'tag' : 'Volume' },
        'Persistent?' : {'type' : 'DP_VARCHAR2CHOICE', 'unit' : None, 'tag' : 'P','choices' : ['Y','N','P*','Uncertain']},
        'Bioaccumulative?' :{'type' : 'DP_YN', 'unit' : None, 'tag' : 'B' },
        'Toxic?' : {'type' : 'DP_YN', 'unit' : None, 'tag' : 'T' } }]
        
        consolidatedDict = {}
        for wdict in wbInfo:
            for key in wdict.keys():
                if key not in consolidatedDict.keys():
                    consolidatedDict[key] = wdict[key]
                else:
                    if consolidatedDict[key] != wdict[key]:
                        print '---'
                        print key,consolidatedDict[key]
                        print key,wdict[key]
                        print '---'
                        exit()
        
        return consolidatedDict

    
    def getXml(self):
        endpointXml = self.makeFloat('ECOSAR Acute fish endpoint (mg/L)',str(self.unit.currentText()),str(self.valueLineEdit.text()))
        endpointXml.tag = 'Endpoint'
        
        endpointNameXml = self.makeQVarchar2('endpointName',str(self.endpointNameCombo.currentText()).encode('ascii'))
        endpointXml.append(endpointNameXml)

        durationXml = self.makeQFloat('duration',str(self.durationUnitCombo.currentText()),str(self.durationLineEdit.text()))
        endpointXml.append(durationXml)

        acuteXml = self.makeQVarchar2('acuteOrChronic',str(self.acuteChronicCombo.currentText()))
        endpointXml.append(acuteXml)

        measurementTypeXml = self.makeQVarchar2('measurementType',str(self.measuredCombo.currentText()))
        endpointXml.append(measurementTypeXml)

        trophicLevelXml = self.makeQVarchar2('trophicLevel',str(self.trophicCombo.currentText()))
        endpointXml.append(trophicLevelXml)

        speciesXml = self.makeQVarchar2('Species',str(self.speciesCombo.currentText()))
        endpointXml.append(speciesXml)
        
        pivotalXml = self.makeQYN('PIVOTAL',str(self.pivotalCombo.currentText()))
        endpointXml.append(pivotalXml)

        activeXml = self.makeQYN('ACTIVE','Y')
        endpointXml.append(activeXml)

        return etree.tostring(endpointXml,pretty_print=True)

class endpointDialog(QtGui.QDialog):
    """ECOSAR dialog."""
    def __init__(self,parent):
        super(endpointDialog,self).__init__(parent)
        self.ew = endpointWidget()    
        self.initUI()

    def initUI(self):
        """Initialise user interface."""
        okButton = QtGui.QPushButton('OK')
        okButton.clicked.connect(self.accept)
        cancelButton = QtGui.QPushButton('Cancel')
        cancelButton.clicked.connect(self.reject)

        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(okButton)
        hbox.addWidget(cancelButton)
 

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ew)
        vbox.addLayout(hbox)

        self.setLayout(vbox)
        self.show()


        
#def main():
#    #create application
#    app = QtGui.QApplication(sys.argv)
#
#    #make instance of example class
#    #ex = ecosarWidget()
#    ed = ecosarDialog(None)
#    edStatus = ed.exec_()
#    if edStatus:
#        print ed.ew.getXml()
#    exit()
#
#
#
#
#    #start event handling (start app)
#    sys.exit(app.exec_())
#
#if __name__== '__main__':
#    main()
#













class ecosarDialogWrapper(reportProvider):
    """Report provider wrapper for the ECOSAR widget."""

    def __init__(self,inputTuple):
        super(ecosarDialogWrapper,self).__init__(inputTuple)
        self.reportType = 'Add endpoint'


    def displayReport(self):
        """Display report."""


        ed = endpointDialog(self.parent)
        edStatus = ed.exec_()
        if edStatus:
            self.chemProxyModel = self.parent.cv.chemProxyModel
            self.dataPointProxyModel = self.parent.cv.dataPointProxyModel
            chemti = self.chemIdx.internalPointer()


            self.model.beginResetModel()
            self.chemProxyModel.beginResetModel()
            self.dataPointProxyModel.beginResetModel()
            
            newXml = ed.ew.getXml()
            print newXml
            newXmlElement = etree.fromstring(newXml)

            self.chemIdx.internalPointer().element.append(newXmlElement)
            self.model.treestep(newXmlElement,self.chemIdx.internalPointer())

            #with open('ddddd.xml','w') as f:
            #    f.write(etree.tostring(self.model.rootItem.element,pretty_print=True))

            self.model.endResetModel()
            self.chemProxyModel.endResetModel()
            self.dataPointProxyModel.endResetModel()

            self.parent.cv.chemView.setRootIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(0,0,self.model.rootItem.child(0))))
            self.parent.cv.chemSelectionModel.setCurrentIndex(self.chemProxyModel.mapFromSource(self.model.createIndex(chemti.row(),0,chemti)),QtGui.QItemSelectionModel.SelectCurrent)
            self.parent.cv.dataPointSelectionModel.setCurrentIndex(self.dataPointProxyModel.mapFromSource(self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1])),QtGui.QItemSelectionModel.SelectCurrent)
            newDataPointIndex = self.model.createIndex(chemti.childItems[-1].row(),0,chemti.childItems[-1])
            self.model.dataChanged.emit(newDataPointIndex,newDataPointIndex)





#    def getDpData(self,dpName,active=False):
#        try:
#            shortName = self.model.structureDict[dpName]['longName']
#
#            xml = self.chemIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'xmlelement')].toPyObject()
#            print './'+shortName
#        except:
#            return []
#
#
#    def getQData(self,qName,active=False):
#        pass
