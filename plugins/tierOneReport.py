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
from reportDialog import reportDialog
import datetime
from lxml import etree
import formatNumbers
import organiseSVGNoRect
import openbabel

class reportTypeOne(reportProvider):
    """Tier I report"""

    def __init__(self,inputTuple):
        super(reportTypeOne,self).__init__(inputTuple)
        self.reportType = 'Tier I Report'

    def makeReport(self):
        """Makes report."""
        import formatNumbers #code to format numbers
        import pickle
        #from time import ctime
        from time import strftime
        from pprint import pprint
        
        def makePanels(panel):
            
            #Dictionary of row names (keys) and HTML for the tooltip (values)
            tooltipDict = {'Persistent?':'<p><b>Persistent if </p> <p>half-life in water is > 60 days </p> or half-life in air > 2 days</p></b></p>', \
                       'Bioaccumulative?':'<p><b>Bioaccumulative if <p>BMF > 1</p><p>or BAF &ge; 2000</p><p> or BCF &ge; 2000 </p><p>or log Kow &ge; 4.2</b></p></p>', \
                       'Toxic?':'<p><b>Toxic if <p>Acute Endpoint &le; 1 mg/L </p><p>or Chronic Endpoint &le; 0.1 mg/L</p></b></p>'
                       } 

           
            
            """function that creates each HTML panel
            
            Takes a list as input with first item as name and then alternating 
            parameter and parameter value."""

            panelName = panel[0] #takes first member of the list (name of the HTML panel)
            #parameters = panel[0::2] 
            #values = panel[1::2]
            rows = ''
        
            for rowList in panel[1]:

                #insert tooltips for the row names (keys) in tooltipDict
                if rowList[0] in tooltipDict.keys():
                        rowList[0] = '<span class="hover-container">' + rowList[0] +'*' + '<div>' + tooltipDict[rowList[0]] + '</div></span>'
                        
                rows = rows + '<tr><td><b>' + rowList[0] + '</b></td>' + '<td>' + rowList[1] + '</td></tr>' + '\n'
        
            table = '<table class="table table-striped">\n<col style="width:33%"/>\n<col style="width:67%"/>\n' + rows + '</table>' + '\n'
            
            HTMLpanel = '<div class="IMAPcontainer" id="'+panelName+'"><div class="panel panel-primary">\n<div class="panel-heading">\n<h1 class="panel-title">\n<b>' + panelName + \
                        '</b></h1></div>' + table + '</div></div>\n\n'
        
            return HTMLpanel
        
        
        
        ########################################################        
        #                                                      #
        #             Begin makeReport proper.                 #
        #                                                      #
        ########################################################

        #reportData = self.getData()
        self.makeNameConversionDicts()
        reportData = self.recursiveGetData(self.chemIdx.internalPointer())
        #pprint(reportData)
        reportData = self.massageData(reportData)
        #pprint(reportData['Assessment Status'])
        #exit()
        #pprint(reportData)

        #pprint(reportData)

        
        self.shortToLong = { self.model.structureDict[key]['longName'] : key for key in self.model.structureDict.keys()}
        #currentTime = ctime()
        currentTime = strftime("%a %d-%b-%Y %I:%M %p",)

        svgString = str(self.chemIdx.data(QtCore.Qt.UserRole).toMap()[QtCore.QString(u'svgString')].toPyObject())


        #make degradant svgs
        #errorSvg = """<?xml version="1.0"?>
        #<svg version="1.1" id="topsvg"
        #xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
        #xmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="200px" height="200px" viewBox="0 0 100 100">
        #<title>OBDepict</title>
        #<rect x="0" y="0" width="100" height="100" fill="white"/>
        #<text text-anchor="middle" font-size="6" fill ="black" font-family="sans-serif"
        #x="50" y="98" ></text>
        #<g transform="translate(0,0)">
        #<svg width="100" height="100" x="0" y="0" viewBox="0 0 80 80"
        #font-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="2"  stroke-linecap="round">
        #<text x="5" y="75" fill="red"  stroke="red" stroke-width="1" font-size="96" >X</text>
        #</svg>
        #</g>
        #</svg>"""

        #degradantSmiles = self.chemIdx.internalPointer().element.xpath('./degradantSMILES/DP_CLOB1/text()')

        #if len(degradantSmiles) >0:
        #    degSvgList = []
        #    obconv = openbabel.OBConversion()
        #    obconv.SetInAndOutFormats('smi','svg')
        #    for smiles in degradantSmiles:
        #        mol = openbabel.OBMol()
        #        readOk=obconv.ReadString(mol,smiles)
        #        if readOk == False:
        #            degSvgList.append(errorSvg)
        #        else:
        #            svg = obconv.WriteString(mol)
        #            degSvgList.append(svg)

        #    p = organiseSVGNoRect.packer(degSvgList)
        #    degradantSvg = p.makeSVG('best',False)
        #    degradantSvg = etree.tostring(degradantSvg,pretty_print=True)
        #else:
        #    degradantSvg = None



        #datalist = getFromExcel.getfromexcel('d:/testfile.xlsm')
        #datalist = pickle.load(open("datalist.p","rb"))
        #
        #rownumber = 20 #row of chemical you wish to make the report
        #chemicalRef = rownumber - 2
        
        #Dictionary of Data Types
        dictDataType = {
        
        'CAS-RN' : 'Integer' ,
        'CAS-RN (hyphenated)' : 'String' ,
        'NICNAS ID' : 'Integer' ,
        'AICS Name' : 'String' ,
        'Common name' : 'String' ,
        'IMAP Status' : 'String' ,
        'Chemical class' : 'String' ,
        'High concern use for environment?' : 'Boolean' ,
        'Perfluorinated? (NICNAS master list 12 Oct 2012)' : 'Boolean' ,
        'OECD HPV' : 'Boolean' ,
        'Montreal' : 'Boolean' ,
        'SGG' : 'Boolean' ,
        'Rotterdam' : 'Boolean' ,
        'Stockholm' : 'Boolean' ,
        'REACH (SVHCs)' : 'Boolean' ,
        'EDC (US EPA)' : 'Boolean' ,
        'EDC (Europe)' : 'Boolean' ,
        'On DSL' : 'Boolean' ,
        'DSL P' : 'Boolean' ,
        'DSL B' : 'Boolean' ,
        'DSL iT' : 'Boolean' ,
        'Persistent?' : 'Boolean' ,
        'Bioaccumulative?' : 'Boolean' ,
        'Toxic?' : 'Boolean' ,
        'Default volume? (8 May 2013)' : 'Boolean' ,
        'PEC(river) (ug/L ; 8 May 2013)' : 'Float' ,
        'Release Mitigation Factor' : 'Float' ,
        'Mitigated PEC(river) (ug/L)' : 'Float' ,
        'Pivotal endpoint (mg/L)' : 'Float' ,
        'Assessment factor' : 'Integer' ,
        'PNEC (ug/L)' : 'Float' ,
        'RQ' : 'Float' ,
        #'Total No. of Uses and Top 3 Uses (SPIN UC62)' : 'String' ,
        'SPIN UC62 - Total Tonnage' : 'Integer' ,
        #'Total No. of Uses - Detailed SPIN' : 'String' ,
        'Comments' : 'String' ,
        'Assessment conclusion' : 'String' ,
        'Assessment Status' : 'String' ,
        'Molecular weight (g/mol)' : 'Float' ,
        'Molecular formula' : 'String' ,
        'SMILES String' : 'String' ,
        'log Kow' : 'Float' ,
        'Water solubility  (mg/L)' : 'Float' ,
        'Melting point (deg C)' : 'Float' ,
        'Vapour pressure (Pa)' : 'Float' ,
        'pKa' : 'Float' ,
        'Ionisable in the environment?' : 'Boolean' ,
        'ECOSAR Acute fish endpoint (mg/L)' : 'Float' ,
        'Fish SAR' : 'String' ,
        'ECOSAR Acute invertebrate endpoint (mg/L)' : 'Float' ,
        'Invertebrate SAR' : 'String' ,
        'ECOSAR Acute algae endpoint (mg/L)' : 'Float' ,
        'Algae SAR' : 'String' ,
        'Batch no.' : 'Integer' ,
        'CAS Valid?' : 'Boolean' ,
        'Fish endpoint (mg/L)' : 'Float' ,
        'Invertebrate endpoint (mg/L)' : 'Float' ,
        'Algae endpoint (mg/L)' : 'Float' ,
        'Pivotal endpoint type' : 'String' ,
        'Primary halflife in water (days; Catalogic 301C)' : 'String' ,
        'BOD (301C)' : 'Float' ,
        'Obs. BOD (301C)' : 'Float' ,
        'SIAR' : 'Boolean' ,
        'Phys Chem Notes' : 'String' ,
        'Reasonable worst case exposure scenario appropriate?' : 'Boolean' ,
        'Reasons for Categorisation (P)' : 'String' ,
        'Reasons for Categorisation (B)' : 'String' ,
        'Assessment Factor Notes' : 'String' ,
        'Internal Notes (not for publication)' : 'String' ,
        'Assessor' : 'String' ,
        'Volume (T) (8 May 2013)' : 'Float' ,
        'Toxicity Notes' : 'String' ,
        'Group' : 'String' ,
        'Fish ECOSAR (mg/L; Neutral Organics SAR)' : 'Float' ,
        'Daphnia ECOSAR (mg/L; Neutral Organics SAR)' : 'Float' ,
        'Algae ECOSAR (mg/L; Neutral Organics SAR)' : 'Float' ,
        'NICNAS Excluded Use' : 'Boolean' ,
        'AHVICL Uses (2006)' : 'String' ,
        'AHVICL Threshold Range (2006)' : 'String'
        
        
        }
        
        
        
        
        #Blueprint to make the HTML panels (in the form of a list of lists). First item is panel name followed by Column Names from Organics Spreadsheet.
        
        panelInfo = {'Key Information': ['RQ','Release Mitigation Factor','Volume (T) (8 May 2013)','Mitigated PEC(river) (ug/L)','PNEC (ug/L)','Persistent?','Bioaccumulative?','Toxic?','Comparison with Canada DSL'],
        'Assessment Summary': ['Assessment conclusion','Comments','Internal Notes (not for publication)','Group','Assessment Status','Batch no.','Assessor'],
        'Use Information' : ['Reasonable worst case exposure scenario appropriate?','AHVICL Uses (2006)','AHVICL Threshold Range (2006)','SPIN UC62 - Total Tonnage','SPIN Uses (UC62)','Detailed SPIN Uses'],
        'Chemical ID Information' : ['NICNAS ID','AICS Name','Common name','IMAP Status','Chemical class','SMILES String','Molecular formula','Molecular weight (g/mol)'],
        'Preassessment Profile' : ['High concern use for environment?','Perfluorinated? (NICNAS master list 12 Oct 2012)','Montreal','SGG','Rotterdam','Stockholm','REACH (SVHCs)','EDC (US EPA)','EDC (Europe)','On DSL','DSL P','DSL B','DSL iT'],
        'Physico-chemical Properties' : ['log Kow','Water solubility  (mg/L)','Melting point (deg C)','Vapour pressure (Pa)','pKa','Ionisable in the environment?','Phys Chem Notes'],
        'PBT Categorisation' :['Persistent?','Obs. BOD (301C)','BOD (301C)','Primary halflife in water (days; Catalogic 301C)','Reasons for Categorisation (P)','Bioaccumulative?','Reasons for Categorisation (B)','Toxic?','Toxicity Notes'],
        'Degradants' : ['Degradant SMILES','Degradants'],
        'ECOSAR - Neutral Organics SAR (mg/L)' :['Fish ECOSAR (mg/L; Neutral Organics SAR)','Daphnia ECOSAR (mg/L; Neutral Organics SAR)','Algae ECOSAR (mg/L; Neutral Organics SAR)'],
        'ECOSAR (mg/L)' : ['ECOSAR Acute fish endpoint (mg/L)','ECOSAR Acute invertebrate endpoint (mg/L)','ECOSAR Acute algae endpoint (mg/L)'],
        'Measured Ecotoxicity Endpoints (mg/L)' : ['Fish endpoint (mg/L)','Invertebrate endpoint (mg/L)','Algae endpoint (mg/L)','Toxicity Notes'],
        'Non-standard Ecotoxicity Endpoints' : ['Other endpoint (mg/L)'],
        'PNEC Calculation' : ['Pivotal endpoint (mg/L)','Pivotal endpoint type','Assessment factor','Assessment Factor Notes','PNEC (ug/L)'],
        'PEC Calculation' : ['Default volume? (8 May 2013)','Volume (T) (8 May 2013)','PEC(river) (ug/L ; 8 May 2013)','Release Mitigation Factor','Mitigated PEC(river) (ug/L)'],
        'Risk Characterisation' : ['RQ','Release Mitigation Factor','Volume (T) (8 May 2013)','Persistent?','Bioaccumulative?','Toxic?']}
        
        panelOrder = ['Key Information','Assessment Summary','Use Information','Chemical ID Information','Preassessment Profile','Physico-chemical Properties','PBT Categorisation','Degradants','ECOSAR - Neutral Organics SAR (mg/L)','ECOSAR (mg/L)','Measured Ecotoxicity Endpoints (mg/L)','Non-standard Ecotoxicity Endpoints', 'PNEC Calculation','PEC Calculation','Risk Characterisation'] 
        
        
        
        
        #DICTIONARY CODE
        #datalistHeader = datalist.pop(0)
        #columnDict = {}
        #
        #for index,item in enumerate(datalistHeader):
        #	columnDict[item] = datalist[chemicalRef][index]
        
        
        
        #print columnDict['pKa']
        #print columnDict.keys()
        #---
        
        #try to find CAS-RN
        
        #Prepare data for HTML
        panelList = []

        dropDownList = ['                    <ul class="dropdown-menu">']

        missing = []
        for key in panelOrder:
            thisPanel = [key,[]]
            #thisPanel = []  
            #thisPanel.append(key)

            #make dropdown
            dropDownList.append('                        <li><a href="#'+key+'">'+key+'</a></li>')

        
            for item in panelInfo[key]:
                #get data
                #data = self.getDpData(item)
                #print data

                #thisPanel[1].append([item,markUp(item,columnDict[item])])
                if item in reportData.keys():
                    print '----'
                    print item,reportData[item]
                    thisPanel[1].append([item,reportData[item]['value']])
                else:
                    thisPanel[1].append([item,''])
                    missing.append(item)
                
                #thisPanel.append(item)
                #thisPanel.append(columnDict[item]) #look up the value of the parameter in the dictionary "columnDict" 
        
            panelList.append(thisPanel)

        dropDownList.append('                    </ul>')
        dropDown = '\n'.join(dropDownList)
        #print dropDown

        
        print 'missing:'
        pprint(missing)
        
        
        
        
        #Components of the HTML
        #Triple double quotes to split string over several lines
        css = """
        <head>
        <meta http-equiv="X-UA-Compatible" content="IE=11;IE=10;IE=9;IE=8; IE=7" />
        <style>

        .wrapper {
        width:600px

        }
        
        .IMAPcontainer
        {
        max-width: 600px;
        min-width: 500px;
        float: left;
        margin-left: 30px;
        margin-bottom: 30px;
        overflow:visible;
        }
        
        .picContainer
        {
        max-width: 300px;
        position: fixed;
        top: 80px;
        right: 50px;
        border: 3px solid gray;
        border-style: outset;
        width: 23%;
        opacity:0.8;
         }
        
        .footer
        {
        clear:both;
        
        margin: 0;
        padding: .5em;
        color: #333;
        background-color: #ddd;
        border-top: 1px solid gray;
        }
         
         table td 
        {
         width:300px;
         overflow:visible;

        -ms-word-break: break-all;
        /*word-break: break-all;*/

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
       
        """
       
        ##The following works if you are connecte to the intertubes.
        #head = """<head> 
        #<link rel="stylesheet" type="text/css" href="http://bootswatch.com/cerulean/bootstrap.min.css">
        #<script type="text/javascript" src="http://code.jquery.com/jquery.min.js"></script>
        #<script type="text/javascript" src="http://bootswatch.com/bower_components/bootstrap/dist/js/bootstrap.min.js"></script>
        #
        #<script type="text/javascript">
        #$(function () {
        #    $("*").tooltip({
        #        placement: "right"
        #		
        #    });
        #});
        #</script>
        #</head>"""


        #Otherwise, you have to read it all in here.
        bsFile = open('plugins/bootstrap.min.css')
        bootstrapCss = bsFile.read()
        bsFile.close()


        bsFile = open('plugins/bootstrap.min.js')
        bootstrapJs = bsFile.read()
        bsFile.close()

        jqFile = open('plugins/jquery.min.js')
        jquery = jqFile.read()
        jqFile.close()
        head = """ 
        <style>
        """+bootstrapCss+ """
        </style>
        <script>
        """+jquery+"""
        </script>
        <script>
        """+bootstrapJs+"""
        </script>
        
        <script type="text/javascript">
        $(function () {
            $("*").tooltip({
                placement: "right"
        		
            });
        });
        </script>
        </head>"""
        
        navbar = """
        
        <div class="nav navbar-inverse navbar-fixed-top">
            <ul class="nav navbar-nav" >
                <li><a href="#">ECAT</a></li>
                <li><a href="#">Tier I Report for """ + reportData['CAS-RN']['value'] +"""</a></li>
        	
                <li class="dropdown">
                    <a href="#" data-toggle="dropdown" class="dropdown-toggle">Report Links<b class="caret"></b></a>
                    """+dropDown+"""
                </li>
                <li><a href="#Assessment Summary">"""+reportData['Assessment Status']['value']+"""</a></li>
               
            </ul>
        </div>        
        
        
        """
        
        footer = """
        
         <div class="nav navbar-inverse navbar-fixed-bottom" style="color:white";>
        Tier I Report <span style="text-align:right";>  --- Created: """ + currentTime + """</span>
        </div>
        
        """
        
        #structure = """
        #
        #<div class="picContainer" >
        #<object width="100%">
        #"""+svgString+"""
        #</object>
        #<p>degradants</p>
        #</div>
        #
        #"""
        
        structure = """
        
        <div class="picContainer" >
        <object width="100%">
        """+svgString+"""
        </object>"""
        
        structure = structure + '\n</div>\n\n'
#        if degradantSvg == None:
#            structure = structure + '\n</div>\n\n'
#        else:
#            structure = structure + '\n<center><b>Degradants</b></center>\n<object width="100%">'+degradantSvg+'\n</object>\n</div>\n'
#
        
        heading = """
        
        <div class="IMAPcontainer">
        <h1>ECAT</h1>
        <h2>Tier I Report for """ + reportData['Display name']['value'] + """</h2>
        <h1><small>CAS-RN: """+ reportData['CAS-RN']['value'] +"""</small></h1>
        </div>
        
        """
        
        space = """
        
        
        <div><h1>&nbsp;</h1></div>
        
        
        """
        
        #write all the components of the HTML code to the file
        panels = ""

        
        
        for item in panelList:
            print 'ITEM: ' + str(item)
            panels = panels + (makePanels(item))

        panels = '<div class="wrapper">\n' + panels + '\n<div>' #the panels are placed into the "wrapper" div so that they align on top of each other (and not side by side when the
                                                                #screen size is changed)
            
        HTMLString = '<html>' + css + '\n' + head + '<body>'+ navbar + space + structure + heading + '\n' + panels + footer + '</body>' + '\n' + '</html>'
        
        return(HTMLString)
        
        
        #print "HTML File: " + fileName + " has been created."



    #def returnReport(self):
    #    return self.reportText

    def recursiveGetData(self,inTi):
        """Recursively gets data from a TreeItem."""
                
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
            print '----'
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

        def addEndpoint(chemData,name,dictData):
            if name not in chemData.keys():
                chemData[name]=[dictData]
            else:
                chemData[name].append(dictData)

        
        #adds HTML markup to selected paramaters.
        def markUp(parameter, value):
            amberValues = ('On DSL', 'EDC (US EPA)','EDC (Europe)')
            redValues = ('Montreal','SGG','Rotterdam','Stockholm','REACH (SVHCs)','Persistent?','Bioaccumulative?','Toxic?','DSL P','DSL B','DSL iT', \
                         'Perfluorinated? (NICNAS master list 12 Oct 2012)', 'High concern use for environment?')
            
            #for uses, splits values according to '#' and adds <p> tag
            #if parameter == 'Total No. of Uses - Detailed SPIN' or parameter =='AHVICL Uses (2006)' or parameter =='Total No. of Uses and Top 3 Uses (SPIN UC62)':
            #        
            #    usesList = value.split("#")
            #    usesList.pop(0) #Pop off the number of uses from Excel spreadsheet. This should be the same as "numberofUses" variable. (i.e. length of usesList).
            #    numberofUses = len(usesList)
        
            #    newList = []
            #    string = ""
        
            #    for i in range(len(usesList)):
            #        usesList[i] = ("<p>" + usesList[i].strip() + "</p>")
            #        string = string + usesList[i]
        
        
            #    string = "<p>Total uses: " + str(numberofUses) + "</p>" + string
            #   
            #    return string
           
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
                if re.match('&lt;',value):
                    #split it
                    splitVal = re.split('&lt;',value.strip())
                    try:
                        if float(splitVal[1]) < 0.01:
                            return "<b><span style=""Color:Green"">&lt;0.01</span></b>"
                        elif float(splitVal[1]) < 1:
                            return "<b><span style=""Color:Green"">" + value + "</span></b>"
                        else:
                            return "<b><span style=""Color:Red"">" + value + "</span></b>"
                    except:
                        return "<b>" + value + "</b>"


                if float(value) < 1.0:
                    if float(value) < 0.01:
                        return "<b><span style=""Color:Green"">&lt;0.01</span></b>"
                    else:
                        return "<b><span style=""Color:Green"">" + formatNumbers.print2dp(value) + "</span></b>"
                else:
                    return "<b><span style=""Color:Red"">" + formatNumbers.print2dp(value) + "</span></b>"
            elif parameter == 'Molecular weight (g/mol)':
                return formatNumbers.print2dp(value)
                return '<b><span style=""Color:#FF7E00"">'+value+'</span></b>'
                    
        
            else:
                return value

        def makeToolTip(key,epDict):

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
                #print 'qKeyList: ' + qKeyList

                qKeyList = sorted(qKeyList, key=getPosition) #sort qualifiers according to the order in "orderList"
                qKeyList.remove('Active')
                toolTipList = []
                for qKey in qKeyList:
                    for item in epDict['children'][qKey]:
                        
##                        if 'http' in item['value']:
##                            
##                            toolTipList.append('<p><b>'+qKey+'</b><a href="'+ item['value']+'">' +item['value']+ '</a></p>') #need to split out text containing "http" and add <a> tags
##                        else:
                            toolTipList.append('<p><b>'+qKey+'</b> '+item['value']+'</p>')
                            #if qKey in shortList:
                            #    shortListQuals += " " + item['value']

                
                #toolTip = ' data-html="true" title="'+'\n'.join(toolTipList)+'"'
                toolTip = '\n'.join(toolTipList)


                #toolTip = ' data-html="true" title="'+epDict['value']+'"'

                #print epDict['children']

                
                
                return toolTip

        def degradantSvgFromSmiles(smiles):


            #make degradant svgs
            errorSvg = """<?xml version="1.0"?>
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

            obconv = openbabel.OBConversion()
            obconv.SetInAndOutFormats('smi','svg')
            mol = openbabel.OBMol()
            readOk=obconv.ReadString(mol,smiles)
            if readOk == False:
                degSvg=errorSvg
            else:
                degSvg = obconv.WriteString(mol)

            p = organiseSVGNoRect.packer([degSvg])
            degradantSvg = p.makeSVG('best',False)
            degradantSvg = etree.tostring(degradantSvg,pretty_print=True)

            return degradantSvg


        ##############################################
        #            Begin massageData               #
        ##############################################

        #print chemData
       
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


        #specialSet = set(['CAS-RN','Common name','Endpoint'])

        endpointSet = set([
            'Fish ECOSAR (mg/L; Neutral Organics SAR)',\
            'Daphnia ECOSAR (mg/L; Neutral Organics SAR)',\
            'Algae ECOSAR (mg/L; Neutral Organics SAR)',\
            'ECOSAR Acute fish endpoint (mg/L)',\
            'ECOSAR Acute invertebrate endpoint (mg/L)',\
            'ECOSAR Acute algae endpoint (mg/L)',\
            'Fish endpoint (mg/L)',\
            'Invertebrate endpoint (mg/L)',\
            'Algae endpoint (mg/L)',\
            'Other endpoint (mg/L)' ])
        specialSet = set(['CAS-RN',\
            'Common name',\
            'Endpoint (mg/L)',\
            'Pivotal endpoint (mg/L)',\
            'Pivotal endpoint type',
            'SPIN UC62 - Total Tonnage',\
            'Comparison with Canada DSL',
            'Degradant SMILES',
            'Degradants']) | endpointSet

        #do CAS-RN
        key = 'CAS-RN'
        if key in chemData.keys():
            if len(chemData[key]) == 1:
                chemData[key]={'value': hyphenate(chemData[key][0]['value'])}
            elif len(chemData[key])>1:
                #join with a comma
                chemData[key]={'value':', '.join([hyphenate(item['value']) for item in chemData[key]])}
        else: # doesn't exist.
            chemData[key]={'value':''}


        #do common name
        key = 'Common name'
        if key in chemData.keys():
            if len(chemData[key]) == 1:
                chemData[key]={'value': chemData[key][0]['value']}
            elif len(chemData[key])>1:
                #join with a comma
                chemData[key]={'value':', '.join([item['value'] for item in chemData[key]])}

        #endpoints
        #pprint(chemData['Endpoint (mg/L)'])
        pivotalList=[]
        pivotalTypeList=[]
        if 'Endpoint (mg/L)' in chemData.keys():
            for ep in chemData['Endpoint (mg/L)']:
                pivotalFlag = False

                #pprint(ep['children'])

                #Categorise. Must be one of
                #Fish ECOSAR (mg/L; Neutral Organics SAR)
                #Daphnia ECOSAR (mg/L; Neutral Organics SAR)
                #Algae ECOSAR (mg/L; Neutral Organics SAR)
                #ECOSAR Acute fish endpoint (mg/L)
                #ECOSAR Acute invertebrate endpoint (mg/L)
                #ECOSAR Acute algae endpoint (mg/L)
                #Fish endpoint (mg/L)
                #Invertebrate endpoint (mg/L)
                #Algae endpoint (mg/L)
                #Other endpoint (mg/L)

                #Is it pivotal? if so, append '(Pivotal endpoint)' to the value.
                #And save some stuff for later use.
                if 'Pivotal endpoint?' in ep['children'].keys():
                    if ep['children']['Pivotal endpoint?'][-1]['value'] == 'Y':
                        pivotalFlag = True
                        pivotalList.append(ep['value'])
                        ep['value'] = ep['value']+ ' (Pivotal Endpoint)'
                else:
                    print 'Should be at least  one pivotal endpoint qualifier.'
                    addEndpoint(chemData,'Other endpoint (mg/L)',ep)

                if 'Test duration' in ep['children'].keys() and 'Endpoint name' in  ep['children'].keys():
                    ep['value'] = ep['value'] + ' (' + ep['children']['Test duration'][-1]['value'] + ' ' + ep['children']['Endpoint name'][-1]['value'] + ')'
                   



                #get endpoint type
                if 'Measurement type (Measured, Calculated, Read-across, Analogue)'\
                        in ep['children'].keys():
                    if len(ep['children']['Measurement type (Measured, Calculated, Read-across, Analogue)']) != 1:
                        print 'Should be exactly one Measurement type qualifier'
                        measurementType = 'Other'
                    else:
                        measurementType = ep['children']['Measurement type (Measured, Calculated, Read-across, Analogue)'][0]['value']
                else:
                    print 'Should be exactly one Measurement type qualifier'
                    measurementType = 'Other'
                
                print measurementType

                if pivotalFlag == True:
                    pivotalTypeList.append(measurementType)


                #ECOSAR should be removed from data (to do).
                if measurementType == 'Calculated' or measurementType == 'ECOSAR':

                    #get SAR qualifier(s), if they exist. Should be only one.
                    if 'SAR' in ep['children'].keys():
                        #It's calculated.
                        if len(ep['children']['SAR']) != 1:
                            print 'Should be exactly one SAR qualifier.'
                            addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                            continue
                    else:
                        print 'Should be exactly one SAR qualifier.'
                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                        continue
                        
                    SAR = ep['children']['SAR'][0]['value']
                    ecosarMatch = re.match('.*ecosar.*',SAR,re.IGNORECASE)
                    if ecosarMatch:
                        #it's an ecosar endpoint.
                                
                        neutralOrganicsMatch = re.match('.*neutral[ \t]+organic.*',SAR,re.IGNORECASE)
                        if neutralOrganicsMatch:
                            #it's a neutral organics ecosar endpoint
                            if 'Trophic Level' in ep['children'].keys():
                                if len(ep['children']['Trophic Level']) != 1:
                                    print 'Should be exactly one Trophic Level qualifier'
                                    addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                                    continue
                            else:
                                print 'Should be exactly one Trophic Level qualifier'
                                addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                                continue

                    
                            trophicLevel = ep['children']['Trophic Level'][0]['value']

                            if trophicLevel == 'Fish':
                                #it's a fish neutral organics ecosar endpoint
                                addEndpoint(chemData,'Fish ECOSAR (mg/L; Neutral Organics SAR)',ep)
                            elif trophicLevel == 'Invertebrate':
                                #it's an invertebrate neutral organics ecosar endpoint
                                addEndpoint(chemData,'Daphnia ECOSAR (mg/L; Neutral Organics SAR)',ep)
                            elif trophicLevel == 'Algae':
                                #it's an algae neutral organics ecosar endpoint
                                addEndpoint(chemData,'Algae ECOSAR (mg/L; Neutral Organics SAR)',ep)
                            else:
                                #it's some whacko calculated ecosar endpoint.
                                addEndpoint(chemData,'Other endpoint (mg/L)',ep)

                        else:
                            #It's an ecosar endpoint, but not neutral organics.
                            
                            #Is it acute or chronic?
                            if 'Acute or chronic endpoint?' in ep['children'].keys():
                                if len(ep['children']['Acute or chronic endpoint?']) not in set([0,1]):
                                    print 'Must be zero or one \'Acute or chronic endpoint?\'.'

                                    addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                                    continuek


                                if ep['children']['Acute or chronic endpoint?'][0]['value'] == 'acute':
                                    #acute ecosar endpoint. Get trophic level.

                                    if 'Trophic Level' in ep['children'].keys():
                                        if len(ep['children']['Trophic Level']) != 1:
                                            print 'Should be exactly one Trophic Level qualifier'
                                            addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                                            continue
                                    else:
                                        print 'Should be exactly one Trophic Level qualifier'
                                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                                        continue

                    
                                    trophicLevel = ep['children']['Trophic Level'][0]['value']


                                    if trophicLevel == 'Fish':
                                        #it's a fish acute ecosar endpoint
                                        addEndpoint(chemData,'ECOSAR Acute fish endpoint (mg/L)',ep)
                                    elif trophicLevel == 'Invertebrate':
                                        #it's an invertebrate acute ecosar endpoint
                                        addEndpoint(chemData,'ECOSAR Acute invertebrate endpoint (mg/L)',ep)
                                    elif trophicLevel == 'Algae':
                                        #it's an algae acute ecosar endpoint
                                        addEndpoint(chemData,'ECOSAR Acute algae endpoint (mg/L)',ep)
                                    else:
                                        #it's some whacko calculated ecosar endpoint.
                                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)

                                else:
                                    #chronic ecosar endpoint.
                                    addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                            else:
                                #No acute or chronic qualifier.
                                addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                    else:
                        #it's a calculated endpoint of unknown origin.
                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                elif measurementType == 'Measured':
                    #it's a measured endpoint.
                   
                    if 'Trophic Level' in ep['children'].keys():
                        if len(ep['children']['Trophic Level']) != 1:
                            print 'Should be exactly one Trophic Level qualifier'
                            addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                            continue
                    else:
                        print 'Should be exactly one Trophic Level qualifier'
                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)
                        continue

                    
                    trophicLevel = ep['children']['Trophic Level'][0]['value']


                    #get trophic level.
                    trophicLevel = ep['children']['Trophic Level'][0]['value']
                    if trophicLevel == 'Fish':
                        #it's a measured fish endpoint
                        addEndpoint(chemData,'Fish endpoint (mg/L)',ep)
                    elif trophicLevel == 'Invertebrate':
                        #it's a measured invertebrate endpoint
                        addEndpoint(chemData,'Invertebrate endpoint (mg/L)',ep)
                    elif trophicLevel == 'Algae':
                        #it's a measured algae endpoint
                        addEndpoint(chemData,'Algae endpoint (mg/L)',ep)
                    else:
                        #it's some other measured endpoint.
                        addEndpoint(chemData,'Other endpoint (mg/L)',ep)

                else:
                    #add to the Other Endpoint section.
                    addEndpoint(chemData,'Other endpoint (mg/L)',ep)

        #Get total spin tonnage.
        if 'SPIN Uses (UC62)' in chemData.keys():
            chemXml = chemData['SPIN Uses (UC62)'][0]['xmlelement'].getparent()
            #volumes = [x.text for x in chemXml.xpath('./SPINUsesUC62/DP_QUANT_VALUE1')]
            volume = sum([float(x.text) for x in chemXml.xpath('./SPINUsesUC62/spinVolume/Q_QUANT_VALUE1')])
            #pprint(chemData['SPIN Uses (UC62)'])
            #exit()
            chemData['SPIN UC62 - Total Tonnage'] = {'value': str(volume) + ' tonne', 'children': {}}

        if 'On DSL' not in chemData.keys():
            chemData['Comparison with Canada DSL'] = {'value': 'Need to check if the chemical is on the Canadian DSL.'}
        else:
            onDslValue = ', '.join([ x['value'] for x in chemData['On DSL']])
            if onDslValue == 'Y':
                #compare DSL categorisations with IMAP PBT categorisations
                if 'DSL P' in chemData.keys():
                    DSLP = ', '.join([ x['value'] for x in chemData['DSL P']])
                else:
                    DSLP = None

                if 'DSL B' in chemData.keys():
                    DSLB = ', '.join([ x['value'] for x in chemData['DSL B']])
                else:
                    DSLB = None


                if 'DSL iT' in chemData.keys():
                    DSLiT = ', '.join([ x['value'] for x in chemData['DSL iT']])
                else:
                    DSLiT = None


                if 'Persistent?' in chemData.keys():
                    IMAPP = ', '.join([ x['value'] for x in chemData['Persistent?']])
                else:
                    IMAPP = None

                if 'Bioaccumulative?' in chemData.keys():
                    IMAPB = ', '.join([ x['value'] for x in chemData['Bioaccumulative?']])
                else:
                    IMAPB = None
                    

                if 'Toxic?' in chemData.keys():
                    IMAPT = ', '.join([ x['value'] for x in chemData['Toxic?']])
                else:
                    IMAPT = None


                if DSLP == IMAPP and DSLB == IMAPB and DSLiT == IMAPT:
                    chemData['Comparison with Canada DSL'] = {'value': 'The domestic PBT categorisation is consistent with the DSL PBiT categorisation.'}

                else:
                    chemData['Comparison with Canada DSL'] = {'value': '<b><span style="Color:#FF7E00">The domestic PBT categorisation is not consistent with the DSL PBiT categorisation.</span></b>'} 
                #print DSLP,DSLB,DSLiT,IMAPP,IMAPB,IMAPT
            elif onDslValue == 'N':
                chemData['Comparison with Canada DSL'] = {'value': 'Not on Canadian DSL.'}
            else:
                chemData['Comparison with Canada DSL'] = {'value': onDslValue+' Need to check manually.'}

        if 'Degradant SMILES' in chemData.keys():
            picList = []
            smilesList = []
            for picDict in chemData['Degradant SMILES']:
                if picDict['children'].keys() == ['Active']:
                    picList.append(degradantSvgFromSmiles(picDict['value'])+'<br/>')
                    smilesList.append(markUp(key,picDict['value'])+'<br/>')
                else:
                    toolTip =makeToolTip(key,picDict)
                    picList.append('<span class="hover-container">'+degradantSvgFromSmiles(picDict['value'])+'<div>'+toolTip+'</div></span><br/>')
                    smilesList.append('<span class="hover-container"><u>'+markUp(key,picDict['value'])+'</u><div>'+toolTip+'</div></span><br/>')
            chemData['Degradants']={'value' : '\n'.join(picList)}
            chemData['Degradant SMILES']={'value' : '\n'.join(smilesList)}
            
            #print chemData['Degradants']
            #print chemData['Degradant SMILES']
            #exit()

        #Make pivotal endpoint and pivotal endpoint type strings for PNEC 
        #calculation.
        chemData['Pivotal endpoint (mg/L)']={}
        chemData['Pivotal endpoint type']={}

        chemData['Pivotal endpoint (mg/L)']['value'] = ', '.join(pivotalList)
        chemData['Pivotal endpoint type']['value'] = ', '.join(pivotalTypeList)

        #make display strings for endpoints
        for key in endpointSet:
            if key in chemData.keys():
                epList = []
                for epDict in chemData[key]:
                    if epDict['children'].keys() == ['Active']:
                        epList.append(markUp(key,epDict['value'])+'<br/>')
                    else:
                        toolTip = makeToolTip(key,epDict)
                        #epList.append('<p'+toolTip+'>'+markUp(key,epDict['value'])+'</p>')
                        epList.append('<span class="hover-container"><u>'+markUp(key,epDict['value'])+'</u><div>'+toolTip+'</div></span><br/>')
                chemData[key]={'value' : '\n'.join(epList)}
        

        #make display strings for everything else

        for key in set(chemData.keys())-specialSet:
            if key in chemData.keys():
                epList = []
                for epDict in chemData[key]:
                    if epDict['children'].keys() == ['Active']:
                        epList.append(markUp(key,epDict['value'])+'<br/>')
                    else:
                        toolTip = makeToolTip(key,epDict)
                        epList.append('<span class="hover-container"><u>'+markUp(key,epDict['value'])+'</u><div>'+toolTip+'</div></span><br/>')
                chemData[key]={'value' : '\n'.join(epList)}


        #Do some final hacks.
        if 'Common name' in chemData.keys():
            chemData['Display name'] = {'value': chemData['Common name']['value']}
        elif 'AICS Name' in chemData.keys():
            chemData['Display name'] = {'value' : chemData['AICS Name']['value']}
        else:
            chemData['Display name'] = {'value' :'-'}

        #add blank assessment status if it isn't already there.
        if 'Assessment Status' not in chemData.keys():
            chemData['Assessment Status'] = {'value': '-'}

        

        return chemData

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


        inputHtml = self.makeReport()
        #outFile = open('out.html','w')
        #outFile.write(inputHtml)
        #outFile.close()
        htmlDialog = reportDialog(self.parent,inputHtml)
        htmlDialog.exec_()

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
