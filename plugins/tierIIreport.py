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
from PyQt4 import QtCore, QtGui,QtWebKit
from pprint import pprint
import re
from reportDialog import textDialog
import datetime
from lxml import etree
import formatNumbers
import pickle
from time import ctime
from pprint import pprint
import makePublicationSvg

class tierII(reportProvider):
    """Tier II report"""

    def __init__(self,inputTuple):
        super(tierII,self).__init__(inputTuple)
        self.reportType = 'Tier II Report'

    def makeReport(self):
        """Makes report."""
        def markUpMF(value):
            subTag = False #True if <sub> tag has been added
            s = ""

            for i in range(len(value)):

                if value[i].isdigit() and subTag == False:
                    s = s + '<sub>' + value[i]
                    subTag = True
                elif not(value[i].isdigit()) and subTag == True:
                   s = s + '</sub>' + value[i]
                   subTag = False
                else:
                   s = s + value[i]

            if subTag == True:
                s = s + '</sub>'

            return s


        def makeTableEntry(ti,key,viewKey):
            if key in ti.keys():
                values = [x['value'] for x in ti[key]]
                if key == 'Molecular formula':
                    for index in range(len(values)):
                        values[index] = markUpMF(values[index])
                elif key == 'SMILES String' and viewKey == 'Structural formula':
                    casrns = [x['value'] for x in ti['CAS-RN']]
                    casrns = list(set(casrns))
                    casrnStr = '_'.join(casrns)+'.jpg'

                    htmlStructure = makePublicationSvg.makeHtml(casrnStr,makePublicationSvg.combinedSvgFromSmiles(values))
                    values = [htmlStructure]

                    

            else:
                values = []



            if len(viewKey) < 27:
                entry = '| '+viewKey+' '*(27-len(viewKey))+'| '
            else:
                entry = '| '+viewKey+'| '
            
            if key in ti.keys():
                entry = entry + '\n|                            | '.join(values)
                
                entry = entry + '\n'
            else:
                entry = entry + '\n'

            return entry
        
        
        ########################################################        
        #                                                      #
        #             Begin makeReport proper.                 #
        #                                                      #
        ########################################################

        #reportData = self.getData()
        self.makeNameConversionDicts()
        #reportData = self.recursiveGetData(self.chemIdx.internalPointer())
        reportData = self.getDataAllChemicals()
        #pprint(reportData)
        #reportData = self.massageData(reportData)
        #pprint(reportData['Assessment Status'])
        #exit()
        #pprint(reportData)

        #pprint(reportData)


        
        self.shortToLong = { self.model.structureDict[key]['longName'] : key for key in self.model.structureDict.keys()}
        currentTime = ctime()

        ##Get display name and assessment status
        #if 'Common name' in reportData.keys():
        #    displayName = ', '.join([cnd['value'] for cnd in reportData['Common name']])
        #elif 'IUPAC Name' in reportData.keys():
        #    displayName = ', '.join([cnd['value'] for cnd in reportData['IUPAC Name']])
        #else:
        #    displayName = '-'

        #if 'CAS-RN' in reportData.keys():
        #    casrn = ', '.join([self.hyphenate(x['value']) for x in reportData['CAS-RN']])
        #else:
        #    casrn = '-'



        #svgString = str(self.chemIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgString')].toPyObject())

        #get all casrns in model
        casrns = []
        for chem in reportData:
            #casrns = casrns + [self.hyphenate(x['value']) for x in chem['CAS-RN']]
            casrns = casrns + [int(x['value']) for x in chem['CAS-RN']]
        casrns = set(casrns)
        casrns = list(casrns)
        casrns = ', '.join([self.hyphenate(str(x)) for x in casrns])

        mdString = """[Insert text here]
===========================================================
Title goes here
-----------------------------------------------------------
CAS Registry Numbers: """ + casrns+'\n-----------------------------------------------------------\n\nPreface\n-------\n'

        with open('plugins/preface.txt','r') as f:
            preface = f.read()
        mdString = mdString+preface+'\n'
        mdString = mdString + """Acronyms & Abbreviations
------------------------
Grouping Rationale
------------------

Chemical Identity
-----------------\n"""

        #make a new table for each ti in model
        for ti in reportData:
            mdString = mdString + '|                            |\n|:---------------------------|:-------------------------------\n'
            mdString = mdString + '| CAS RN                     | '+', '.join([self.hyphenate(x['value']) for x in ti['CAS-RN']])+'\n'

            #mdString = mdString + '|Chemical Name              |'+ '\n|                           |'.join(x['value'] for x in ti['AICS Name'])
            mdString = mdString + makeTableEntry(ti,'AICS Name','Chemical Name')
            mdString = mdString + makeTableEntry(ti,'Common name','Synonyms')
            mdString = mdString + makeTableEntry(ti,'SMILES String','Structural formula')
            mdString = mdString + makeTableEntry(ti,'Molecular formula','Molecular formula')
            mdString = mdString + makeTableEntry(ti,'Molecular weight (g/mol)','Molecular weight (g/mol)')
            mdString = mdString + makeTableEntry(ti,'SMILES String','SMILES')
            mdString = mdString + '\n\n<br/>'


        mdString = mdString + """Physical and Chemical Properties
--------------------------------
Import, Manufacture and Use
---------------------------
###Australia

###International

Environmental Regulatory Status
-------------------------------
###Australia

###United Nations

###OECD

###Canada

###European Union

###United States of America


Environmental Exposure
----------------------

Environmental Effects
---------------------

Categorisation of Environmental Hazard
--------------------------------------

Risk Characterisation
---------------------

Key Findings
------------

Recommendations
---------------

"""
        #print ' '.join(x['value'] for x in reportData['Common name'])
        #print ' '.join(x['value'] for x in reportData['Molecular weight (g/mol)'])
        #print ' '.join(x['value'] for x in reportData['Molecular formula'])
        #print ' '.join(x['value'] for x in reportData['SMILES String'])
        #print self.chemIdx.internalPointer().parent().childItems
        


        return(mdString)
    
    def getDataAllChemicals(self):
        def sortKey(ti):
            a =  [int(x['value']) for x in  ti['CAS-RN']]
            a.sort()
            a = tuple(a)
            return a

        chemProxyModel = self.parent.cv.chemProxyModel
        proxyRootIdx = chemProxyModel.mapFromSource(self.chemIdx).parent()
        reportData = []
        for row in xrange(chemProxyModel.rowCount(proxyRootIdx)):
            sourceIdx = chemProxyModel.mapToSource(proxyRootIdx.child(row,0))
            ti = sourceIdx.internalPointer()
            data = self.recursiveGetData(ti)
            data = self.massageData(data)
            reportData.append(data)


        #reportData.sort(key=lambda item: ', '.join([x['value'] for x in  item['CAS-RN']]) )
        reportData.sort(key=sortKey )
        for item in reportData:
            print ', '.join([x['value'] for x in  item['CAS-RN']])

        return reportData



        
    def makeEndpointPanel(self,reportData):
        """Make effects panel.
        
        Three columns, species and organism type on the left, acute/chronic in
        the middle, value on the right.."""

        panelName = 'Effects'

        if 'Endpoint (mg/L)' in reportData.keys():
            epDictList = reportData['Endpoint (mg/L)']
        else:
            epDictList = []
        rows = []
        for epDict in epDictList:
            if epDict['confidential'] == False:
                row = []
                #get species
                if 'Species' in epDict['children'].keys():
                    #Schema says there can be only one.
                    species = '<br/>(<i>'+epDict['children']['Species'][0]['value']+'</i>)'

                else:
                    species = ''

                #get organism type
                if 'Trophic Level' in epDict['children'].keys():
                    organismType = epDict['children']['Trophic Level'][0]['value']
                else:
                    organismType = ''

                row.append(organismType+' '+species)

                #get acute/chronic
                if 'Acute or chronic endpoint?' in epDict['children'].keys():
                    acuteOrChronic = epDict['children']\
                            ['Acute or chronic endpoint?'][0]['value']
                else:
                    acuteOrChronic = ''

                row.append(acuteOrChronic)

                #get test duration
                if 'Test duration' in epDict['children'].keys():
                    testDuration = epDict['children']\
                            ['Test duration'][0]['value']
                else:
                    testDuration = ''

                if 'Endpoint name' in epDict['children'].keys():
                    endpointName = epDict['children']\
                            ['Endpoint name'][0]['value']
                else:
                    endpointName = ''

                if epDict['children'].keys() == ['Active']:
                    value = self.markUp(rowName,epDict['value'])
                else:
                    toolTip = self.makeToolTip('Endpoint (mg/L)',epDict)
                    value = '<span class="hover-container"><u>'\
                            +self.markUp('Endpoint (mg/L)',epDict['value'])+'</u><div>'\
                            +toolTip+'</div></span>'

                row.append(testDuration+' '+endpointName+' = '+value)
                rows.append('<tr><td>'+'</td><td>'.join(row)+'</td></tr>')
            elif epDict['confidential'] == True:
                rows.append('<tr><td style="color:red">CONFIDENTIAL</td><td/><td/></tr>')
            else:
                rows.append('<tr><td/></td><td/><td/></tr>')

        table = '<table class="table table-striped">\n<col style="width:33%"/>\n<col style="width:33%"/>\n<col style="width:33%"/>\n' + '\n'.join(rows) + '</table>' + '\n'
        
        HTMLpanel = '<div class="IMAPcontainer" id="'+panelName+'"><div class="panel panel-primary">\n<div class="panel-heading">\n<h1 class="panel-title">\n<b>' + panelName + \
                    '</b></h1></div>' + table + '</div></div>\n\n'
    
        return HTMLpanel
        
    def makeFatePanel(self,reportData):
        """Make fate panel.
        
        Three columns, species and organism type on the left, acute/chronic in
        the middle, value on the right.."""

        panelName = 'Fate'

        rowName = 'Dissipation time'
        if 'Dissipation time' in reportData.keys():
            epDictList = reportData['Dissipation time']
        else:
            epDictList = []
        rows = []
        for epDict in epDictList:
            if epDict['confidential'] == False:
                row = []
                #get type
                if 'Type' in epDict['children'].keys():
                    #Schema says there can be only one.
                    Type = epDict['children']['Type'][0]['value']

                else:
                    Type = ''
                row.append(Type)

                if epDict['children'].keys() == ['Active']:
                    value = self.markUp(rowName,epDict['value'])
                else:
                    toolTip = self.makeToolTip('Dissipation time',epDict)
                    
                    value = '<span class="hover-container"><u>'\
                            +self.markUp('Dissipation time',epDict['value'])+'</u><div>'\
                            +toolTip+'</div></span>'

                #check to see whether it's 'Not stable to degradation'
                if 'DP_QUAL_VALUE1' in [x.tag for x in epDict['xmlelement']]:
                    row.append(value)
                else:
                    #get name of endpoint (DT50, DT90, etc.)
                    dtName = epDict['xmlelement'].xpath('./Name/Q_QUAL_VALUE1/text()')
                    if len(dtName) ==0:
                        row.append(value)
                    elif len(dtName) == 1:
                        row.append(dtName[0]+' = '+value)
                    else:
                        row.append(value)

                #get comment
                if 'Comment' in epDict['children'].keys():
                    #comment = epDict['children']\
                    #        ['Comment'][0]['value']
                    comment = '<br/>\n'.join(x['value'] for x in epDict['children']\
                            ['Comment'])
                else:
                    comment = ''

                row.append(comment)

                rows.append('<tr><td>'+'</td><td>'.join(row)+'</td></tr>')
            elif epDict['confidential'] == True:
                rows.append('<tr><td style="color:red">CONFIDENTIAL</td><td/><td/></tr>')
            else:
                rows.append('<tr><td/></td><td/><td/></tr>')


        rowList = ['Bioconcentration factor','Soil organic carbon - water distribution coefficient (Koc)','Soil/Water distribution coefficient (Kd)']

        for rowName in rowList:
            if rowName in reportData.keys():
                epDictList = reportData[rowName]
            else:
                epDictList = []
            
            #rows = []
            for epDict in epDictList:
                row = []
                #get type
                #if 'Type' in epDict['children'].keys():
                #    #Schema says there can be only one.
                #    Type = epDict['children']['Type'][0]['value']

                #else:
                #    Type = ''
                row.append(rowName)

                if epDict['children'].keys() == ['Active']:
                    value = self.markUp(rowName,epDict['value'])
                else:
                    toolTip = self.makeToolTip('Dissipation time',epDict)
                    value = '<span class="hover-container"><u>'\
                            +self.markUp('Dissipation time',epDict['value'])+'</u><div>'\
                            +toolTip+'</div></span>'

                row.append(value)

                #get comment
                if 'Comment' in epDict['children'].keys():
                    comment = epDict['children']\
                            ['Comment'][0]['value']
                else:
                    comment = ''

                row.append(comment)

                rows.append('<tr><td>'+'</td><td>'.join(row)+'</td></tr>')


        table = '<table class="table table-striped">\n<col style="width:33%"/>\n<col style="width:33%"/>\n<col style="width:33%"/>\n' + '\n'.join(rows) + '</table>' + '\n'
        
        HTMLpanel = '<div class="IMAPcontainer" id="'+panelName+'"><div class="panel panel-primary">\n<div class="panel-heading">\n<h1 class="panel-title">\n<b>' + panelName + \
                    '</b></h1></div>' + table + '</div></div>\n\n'
    
        return HTMLpanel


    def markUp(self,parameter, value):
        amberValues = ('On DSL', 'EDC (US EPA)','EDC (Europe)')
        redValues = ('Montreal','SGG','Rotterdam','Stockholm','REACH (SVHCs)','Persistent?','Bioaccumulative?','Toxic?','DSL P','DSL B','DSL iT', \
                     'Perfluorinated? (NICNAS master list 12 Oct 2012)', 'High concern use for environment?')
        
        #adds subscripts to molecular formula
        if parameter == 'Molecular formula':
                
            subTag = False #True if <sub> tag has been added
            s = ""
    
            for i in range(len(value)):
    
                if value[i].isdigit() and subTag == False:
                    s = s + '<sub>' + value[i]
                    subTag = True
                elif not(value[i].isdigit()) and subTag == True:
                   s = s + '</sub>' + value[i]
                   subTag = False
                else:
                   s = s + value[i]
    
            if subTag == True:
                s = s + '</sub>'
    
            return s
    
        
        elif (value == 'Yes' or value == 'Y') and parameter in amberValues:#bold and colour amber
                
            return "<b><span style=""Color:#FF7E00"">" + value + "</span></b>"
    
        elif (value == 'Yes' or value == 'Y' or value == 'P*') and parameter in redValues:#bold and colour red
                
            return "<b><span style=""Color:Red"">" + value + "</span></b>"
    
        elif parameter == 'RQ':
            if float(value) < 1.0:
                return "<b><span style=""Color:Green"">" + formatNumbers.print2dp(value) + "</span></b>"
            else:
                return "<b><span style=""Color:Red"">" + formatNumbers.print2dp(value) + "</span></b>"
        elif parameter == 'Molecular weight (g/mol)':
            return formatNumbers.print2dp(value.replace('g/mol','').strip())+' g/mol'
            #return '<b><span style=""Color:#FF7E00"">'+value+'</span></b>'

        elif parameter == 'CAS-RN':
            return self.hyphenate(value)
                
    
        else:
            return value


    def makeToolTip(self,key,epDict):

        orderList = ['Trophic Level','Species','Acute or chronic endpoint?','Endpoint name', \
                     'Test duration', 'Pivotal endpoint?', 'Measurement type (Measured, Calculated, Read-across, Analogue)', 'Comment', 'Citation']

        #shortList = ['Test duration','Endpoint name']
        #shortListQuals = ""

        def getPosition(param):
            if param in orderList:
                return orderList.index(param)
            else:
                return len(orderList) + 1
        
        #print 'makeToolTip: ',key,[epDict]
        if epDict['children'].keys() == ['Active']:
            return ''
        else:
            qKeyList = epDict['children'].keys()

            #sort qualifiers according to the order in "orderList"
            qKeyList = sorted(qKeyList, key=getPosition)
            qKeyList.remove('Active')
            toolTipList = []
            for qKey in qKeyList:
                for item in epDict['children'][qKey]:
                    
                      if 'http' in item['value']:
                          
                          toolTipList.append('<p><b>'+qKey+'</b><a href="'+ item['value']+'">' +item['value']+ '</a></p>') #need to split out text containing "http" and add <a> tags
                      else:
                        toolTipList.append('<p><b>'+qKey+'</b> '+item['value']+'</p>')
            toolTip = '\n'.join(toolTipList)

            return toolTip

    def recursiveGetData(self,inTi):
        """Recursively gets data from a TreeItem."""
 
        def isConfidentialDp(dp):
            try:
                conf1 = dp.element.xpath('./DP_CONFIDENTIAL1_YN/text()')[0]
                #print 'conf1:',conf1
                if conf1 == 'Y':
                    conf1 = True
                elif conf1 == 'N':
                    conf1 = False
                #print 'conf1:',conf1
            except:
                conf1 = None
                    
            try:
                conf2 = dp.element.xpath('./DP_CONFIDENTIAL2_YN/text()')[0]
                #print 'conf2:',conf2
                if conf2 == 'Y':
                    conf2 = True
                elif conf2 == 'N':
                    conf2 = False
                #print 'conf2:',conf2
            except:
                conf2 = None

            if conf1 == None:
                return conf2

            if conf2 == None:
                return conf1

            return conf1 or conf2

        def isConfidentialQ(dp):
            try:
                conf1 = dp.element.xpath('./Q_CONFIDENTIAL1_YN/text()')[0]
                #print 'conf1:',conf1
                if conf1 == 'Y':
                    conf1 = True
                elif conf1 == 'N':
                    conf1 = False
                #print 'conf1:',conf1
            except:
                conf1 = None
                    
            try:
                conf2 = dp.element.xpath('./Q_CONFIDENTIAL2_YN/text()')[0]
                #print 'conf2:',conf2
                if conf2 == 'Y':
                    conf2 = True
                elif conf2 == 'N':
                    conf2 = False
                #print 'conf2:',conf2
            except:
                conf2 = None
                #print 'conf2:',conf2
            
            if conf1 == None:
                return conf2

            if conf2 == None:
                return conf1

            #print conf1,conf2
            return conf1 or conf2

        def tiType(ti):
            chemId = ti['xmlelement'].xpath('./CHEMICAL_ID')
            dpId = ti['xmlelement'].xpath('./DATA_POINT_ID')
            qId = ti['xmlelement'].xpath('./QUALIFIER_ID')

            if len(chemId) == 1 and len(dpId) == 1:
                return 'dp'
            elif len(dpId) == 1 and len(qId) == 1:
                return 'q'
            else:
                #other; either chemical or root.
                return 'o'

        
        #Begin recursiveGetData proper.
                
        contents = {}
        for ti in inTi.childItems:
            tempIdx = self.model.createIndex(ti.row(),0,ti)
            value = str(tempIdx.data(QtCore.Qt.DisplayRole).toString())
            #strip off the bold title used in the viewer.
            value = re.sub('^<b>.*?</b>[ \t]*','',value)

            tag = self.shortToLong[str(ti.element.tag)]
            xml = ti.element

            spinFlag = False
            spinVolume = 0.0
            
            tiDict = {'treeItem': ti, 'xmlelement': xml,'value':value, 'children': self.recursiveGetData(ti)}

            tt = tiType(tiDict)
            #print tt
            if tt == 'dp':
                tiDict['confidential'] = isConfidentialDp(ti)
                #print tiDict['confidential']
            elif tt == 'q':
                tiDict['confidential'] = isConfidentialQ(ti)
                #print tiDict['confidential']



            if tag not in contents.keys():
                contents[tag] = [tiDict]
            else:
                contents[tag].append(tiDict)



        return contents


    def massageData(self,chemData):
        """Massages data gathered from recursiveGetData."""
        def hyphenate(CASRN):
            return CASRN[:-3]+'-'+CASRN[-3:-1]+'-'+CASRN[-1]

        def isActive(dp):
            #find all active qualifiers
            activeQualifiers = dp['children']['Active']


            #are there any 'AUTO' qualifiers? If so, pick the
            #last one in the list.
            mostRecentQualifier = None
            #print '----'
            for item in activeQualifiers:

                if mostRecentQualifier == None:
                    mostRecentQualifier = item
                    continue
                

                #get created date
                createdTimeDate =item['xmlelement'].xpath('./CREATED_DATE')[0].text

                if createdTimeDate == 'AUTO':
                    mostRecentQualifier = item
                    continue
                else:
                    #should be iso 8601 format.
                    if mostRecentQualifier == 'AUTO':
                        continue
                    else:
                        #mostRecentQualifier should be iso 8601 format
                        test = datetime.datetime.strptime(createdTimeDate,"%Y-%m-%dT%H:%M:%S")
                        currentStr = mostRecentQualifier['xmlelement'].xpath('./CREATED_DATE')[0].text
                        current = datetime.datetime.strptime(currentStr,"%Y-%m-%dT%H:%M:%S")
                        if test >= current:
                            mostRecentQualifier = item



            if mostRecentQualifier['value'] == 'Y':
                return True
            elif mostRecentQualifier['value'] == 'N':
                return False
            else:
                print 'massageData: ???'
                exit()



        ##############################################
        #            Begin massageData               #
        ##############################################

        #firstly, get rid of anything that isn't active
        for key in chemData.keys():
            #go from end of list to beginning.
            for i in reversed(range(len(chemData[key]))):
                if isActive(chemData[key][i]) == False:
                    if len(chemData[key]) > 1:
                        del chemData[key][i]
                    elif len(chemData[key]) == 1:
                        del chemData[key]
                        break
                    #pprint(chemData[key])


        return chemData

    
    def hyphenate(self,CASRN):
        return CASRN[:-3]+'-'+CASRN[-3:-1]+'-'+CASRN[-1]

    def makeNameConversionDicts(self):
        """Make dictionaries to convert long to short names and vice-versa."""
        #do datapoint names
        self.shortToLong = { key:self.model.structureDict[key]['longName'] for key in self.model.structureDict.keys()}

        #for completeness
        self.shortToLong['CHEMICAL']='CHEMICAL'
        self.shortToLong['ROOT']='ROOT'

        #do qualifiers
        notOk = []
        for dpKey in self.model.structureDict.keys():
            for qKey in self.model.structureDict[dpKey]['qualifiers'].keys():
                if qKey not in self.shortToLong.keys():
                    self.shortToLong[qKey] = self.model.structureDict[dpKey]['qualifiers'][qKey]['longName']
                else:
                    notOk.append(qKey)

        self.longToShort = { self.shortToLong[key]: key for key in self.shortToLong.keys()}


    def displayReport(self):
        """Display report."""


        inputText = self.makeReport()
        #outFile = open('out.html','w')
        #outFile.write(inputHtml)
        #outFile.close()
        td = textDialog(self.parent,inputText)
        td.exec_()

