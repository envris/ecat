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


from lxml import etree
import sys
import cx_Oracle
from pprint import pprint
import re
import getpass
import datetime
import sqlite3
from PyQt4 import QtGui,QtCore



class uploadfromxml:
    def __init__(self,schema,parentWidget):
        self.parent = parentWidget
        #The db implementation will determine how the connection is set up.
        self.con = cx_Oracle.connect('schema','password','db')
        self.cur = self.con.cursor()
        self.cur.arraysize = 256
        self.allChemicalIDs = []
        self.schema = schema
        
        
    def upload(self,xml):

        #Get Windows user name to put into xml
        author = getpass.getuser()

        if not self.schema.validate(xml):
            print 'xml does not validate!'
            print xml.tag,xml.text
            return

        else:
            print 'xml validates'

            print 'checking VARHCAR2 lengths ...'

            xpathStrings = ['//DP_QUAL_VALUE1/text()','//DP_QUAL_VALUE2/text()',\
                            '//Q_QUAL_VALUE1/text()','//Q_QUAL_VALUE2/text()']
            for xpathString in xpathStrings:
                qualValStrings = xml.xpath(xpathString)
                for item in qualValStrings:
                    if len(item) > 4000:
                        printString = '\"'+str(item) + '\" is longer than 4000 characters.'
                        reply = QtGui.QMessageBox.critical(self.parent,'Danger, Will Robinson!',
                            printString,QtGui.QMessageBox.Ok,QtGui.QMessageBox.Ok)
                        return
                

            #iterate through chemicals
            chemno=0
            for chemical in xml:
                chemno += 1
                
                #check CHEMICAL_ID
                chemical_id = chemical.xpath('./CHEMICAL_ID')
                if len(chemical_id) !=  1:
                    print 'Not exactly one chemical_id!'
                    sys.exit()
                else:
                    chemical_id = chemical_id[0].text
                
                    
                print chemno,'chemical_id: ',chemical_id
                
                if chemical_id == 'AUTO':
                    #upload it, and get new chemical_id
##                    created_by = chemical.xpath('./CREATED_BY')
##                    if len(created_by) != 1:
##                        print 'Not exactly one created_by!'
##                        sys.exit()
##                    else:
##                        created_by = created_by[0].text

                    
                    chemdict = {'CREATED_BY': author}
                    self.cur.execute('insert into CHEMICAL (CHEMICAL_ID,CREATED_BY,CREATED_DATE) '\
                            'values (chemical_id_seq.nextval,:CREATED_BY,SYSDATE)',chemdict)

                    #get CHEMICAL_ID of new entry
                    self.cur.execute('select chemical_id_seq.currval from dual')
                    chemical_id = self.cur.fetchall()
                    if len(chemical_id) != 1:
                        print 'Not exactly one chemical_id!!'
                        sys.exit()
                    else:
                        chemical_id = str(chemical_id[0][0])
##                                      Note that this stuff is unnecessary - this is already in the schema.
##                                else:
##                                        #check to see that chemical_id is an integer
##                                        if re.match('^[0-9]+$',chemical_id) == None:
##                                                print 'chemical_id is not AUTO or an integer.'
##                                                sys.exit()
##                                        #note: if chemical_id != 'AUTO', we already have it stored.
                    
                #now iterate through DATA_POINTs (everything with a child)
                datapoints = chemical.xpath('./*[*]')

                for datapoint in datapoints:
                    #check DATA_POINT_ID
                    data_point_id = datapoint.xpath('./DATA_POINT_ID')
                    if len(data_point_id) != 1:
                        print 'Not exactly one data_point_id!'
                        sys.exit()
                    else:
                        data_point_id = data_point_id[0].text

                    #Note that some  of these (e.g CREATED_DATE) are not used.

                    if data_point_id == 'AUTO':
                        #Note that some  of these (e.g CREATED_DATE) are not used.

                        datapointdict = {'CHEMICAL_ID' : '',
                             'DP_L0_DOMAIN_CODE' : '',
                             'DP_L1_KINGDOM_CODE' : '',
                             'DP_QUAL_VALUE1' : '',
                             'DP_QUANT_VALUE1' : '',
                             'DP_INTEGER1' : '',
                             'DP_UNIT_CODE1' : '',
                             'DP_CONFIDENTIAL1_YN' : '',
                             'DP_CLOB1' : '',
                             'DP_QUAL_VALUE2' : '',
                             'DP_QUANT_VALUE2' : '',
                             'DP_INTEGER2' : '',
                             'DP_UNIT_CODE2' : '',
                             'DP_CONFIDENTIAL2_YN' : '',
                             'CREATED_BY' : ''}
                        datapointdict['DP_L0_DOMAIN_CODE'] = datapoint.tag
                        for item in datapoint.xpath('./*[not(*)]'):
                            if item.text != 'AUTO' and item.tag in datapointdict.keys():
                                datapointdict[item.tag] = item.text
                            elif item.text != 'AUTO' and item.tag not in datapointdict.keys():
                                print 'not in datapointdict: ', item.tag,item.text


    

                        #Set CHEMICAL_ID
                        datapointdict['CHEMICAL_ID']=chemical_id
                        datapointdict['CREATED_BY']=author

                        ##Fix up single quotes in names (need to double up the quotes for Oracle).
                        #datapointdict['DP_CLOB1'] = datapointdict['DP_CLOB1'].replace('\'','\'\'')
                        #pprint(datapointdict)
                        self.cur.execute('insert into DATA_POINT (DATA_POINT_ID, CHEMICAL_ID, DP_L0_DOMAIN_CODE, '\
                                 'DP_L1_KINGDOM_CODE, DP_QUAL_VALUE1, DP_QUANT_VALUE1, DP_INTEGER1, '\
                                 'DP_UNIT_CODE1, DP_CONFIDENTIAL1_YN, DP_CLOB1, DP_QUAL_VALUE2, '\
                                 'DP_QUANT_VALUE2, DP_INTEGER2, DP_UNIT_CODE2,DP_CONFIDENTIAL2_YN, '\
                                 'CREATED_BY, CREATED_DATE) '\
                                 'values (data_point_id_seq.nextval,:CHEMICAL_ID,:DP_L0_DOMAIN_CODE,'\
                                 ':DP_L1_KINGDOM_CODE,:DP_QUAL_VALUE1,:DP_QUANT_VALUE1,:DP_INTEGER1,'\
                                 ':DP_UNIT_CODE1,:DP_CONFIDENTIAL1_YN,:DP_CLOB1,:DP_QUAL_VALUE2,'\
                                 ':DP_QUANT_VALUE2,:DP_INTEGER2,:DP_UNIT_CODE2,:DP_CONFIDENTIAL2_YN,'\
                                 ':CREATED_BY,SYSDATE)',datapointdict)

                        #get new DATA_POINT_ID of new entry
                        self.cur.execute('select data_point_id_seq.currval from dual')
                        data_point_id = self.cur.fetchall()
                        if len(data_point_id) != 1:
                            print 'Not exactly one data_point_id !!'
                            sys.exit()
                        else:
                            data_point_id = str(data_point_id[0][0])

                    #Now iterate through QUALIFIERs (everything with a child)
                    qualifiers = datapoint.xpath('./*[*]')
                    for qualifier in qualifiers:
                        qualifier_id = qualifier.xpath('./QUALIFIER_ID')
                        
                        if len(qualifier_id) != 1:
                            print 'Not exactly one qualifier_id!'
                            sys.exit()
                        else:
                            qualifier_id = qualifier_id[0].text

                        if qualifier_id == 'AUTO':
                            qualdict = {'DATA_POINT_ID' : '',
                                 'Q_L0_DOMAIN_CODE' : '',
                                 'Q_L1_KINGDOM_CODE' : '',
                                 'Q_QUAL_VALUE1' : '',
                                 'Q_QUANT_VALUE1' : '',
                                 'Q_INTEGER1' : '',
                                 'Q_UNIT_CODE1' : '',
                                 'Q_CONFIDENTIAL1_YN' : '',
                                 'Q_CLOB1' : '',
                                 'Q_QUAL_VALUE2' : '',
                                 'Q_QUANT_VALUE2' : '',
                                 'Q_INTEGER2' : '',
                                 'Q_UNIT_CODE2' : '',
                                 'Q_CONFIDENTIAL2_YN' : '',
                                 'CREATED_BY' : ''}
                            qualdict['Q_L0_DOMAIN_CODE'] = qualifier.tag
                            for item in qualifier.xpath('./*[not(*)]'):
                                if item.text != 'AUTO' and item.tag in qualdict.keys():
                                    qualdict[item.tag] = item.text
                                elif item.text != 'AUTO' and item.tag not in qualdict.keys():
                                    print 'not in qualdict: ', item.tag,item.text

                            #set DATA_POINT_ID
                            qualdict['DATA_POINT_ID'] = data_point_id
                            qualdict['CREATED_BY']=author
                            #pprint(qualdict)
                            self.cur.execute('insert into QUALIFIER (QUALIFIER_ID, DATA_POINT_ID, Q_L0_DOMAIN_CODE, '\
                                     'Q_L1_KINGDOM_CODE, Q_QUAL_VALUE1, Q_QUANT_VALUE1, Q_INTEGER1, Q_UNIT_CODE1, '\
                                     'Q_CONFIDENTIAL1_YN, Q_CLOB1, Q_QUAL_VALUE2, Q_QUANT_VALUE2, Q_INTEGER2, '\
                                     'Q_UNIT_CODE2, Q_CONFIDENTIAL2_YN,CREATED_BY, CREATED_DATE) '\
                                     'values (qualifier_id_seq.nextval,:DATA_POINT_ID,:Q_L0_DOMAIN_CODE,'\
                                     ':Q_L1_KINGDOM_CODE,:Q_QUAL_VALUE1,:Q_QUANT_VALUE1,:Q_INTEGER1,:Q_UNIT_CODE1,'\
                                     ':Q_CONFIDENTIAL1_YN,:Q_CLOB1,:Q_QUAL_VALUE2,:Q_QUANT_VALUE2,:Q_INTEGER2,'\
                                     ':Q_UNIT_CODE2,:Q_CONFIDENTIAL2_YN,:CREATED_BY,SYSDATE)',qualdict)
        self.con.commit()
        print 'Upload complete!'
		




class downloadToXml:
    def __init__(self):
        
        self.localDb =  sqlite3.Connection(':memory:',detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.localDb.execute('create table chemical(CHEMICAL_ID INTEGER,CREATED_BY TEXT,CREATED_DATE TIMESTAMP)')
        self.localDb.execute('create table data_point(DATA_POINT_ID INTEGER, CHEMICAL_ID INTEGER, DP_L0_DOMAIN_CODE TEXT, '\
                'DP_L1_KINGDOM_CODE TEXT, DP_QUAL_VALUE1 TEXT, DP_QUANT_VALUE1 REAL, DP_INTEGER1 INTEGER, '\
                'DP_UNIT_CODE1 TEXT, DP_CONFIDENTIAL1_YN TEXT, DP_CLOB1 TEXT, DP_QUAL_VALUE2 TEXT, '\
                'DP_QUANT_VALUE2 REAL, DP_INTEGER2 INTEGER, DP_UNIT_CODE2 TEXT,DP_CONFIDENTIAL2_YN TEXT, '\
                'CREATED_BY TEXT, CREATED_DATE TIMESTAMP)')
        self.localDb.execute('create table qualifier(QUALIFIER_ID INTEGER, DATA_POINT_ID INTEGER, Q_L0_DOMAIN_CODE TEXT, '\
                'Q_L1_KINGDOM_CODE TEXT, Q_QUAL_VALUE1 TEXT, Q_QUANT_VALUE1 REAL, Q_INTEGER1 INTEGER, Q_UNIT_CODE1 TEXT, '\
                'Q_CONFIDENTIAL1_YN TEXT, Q_CLOB1 TEXT, Q_QUAL_VALUE2 TEXT, Q_QUANT_VALUE2 REAL, Q_INTEGER2 INTEGER, '\
                'Q_UNIT_CODE2 TEXT, Q_CONFIDENTIAL2_YN TEXT,CREATED_BY TEXT, CREATED_DATE TIMESTAMP)')

        #The db implementation will determine how the connection is set up.
        self.oracleCon = cx_Oracle.connect('schema','passwd','db')
        oracleCur = self.oracleCon.cursor()

        oracleCur.execute('select chemical_id from chemical')
        self.chemids = oracleCur.fetchall()
        self.chemids = [x[0] for x in self.chemids]

        #populate chemical table
        oracleCur.execute('select chemical_id,created_by,created_date from chemical')

        chemical = []
        for line in oracleCur:
            chemical.append(line)
        #pprint(chemical)

        self.localDb.executemany('insert into chemical values (?,?,?)',chemical)
        del chemical
        print 'chemical table done'

        #populate data_point table
        oracleCur.execute('select DATA_POINT_ID, CHEMICAL_ID, DP_L0_DOMAIN_CODE, '\
                'DP_L1_KINGDOM_CODE, DP_QUAL_VALUE1, DP_QUANT_VALUE1, DP_INTEGER1, '\
                'DP_UNIT_CODE1, DP_CONFIDENTIAL1_YN, DP_CLOB1, DP_QUAL_VALUE2, '\
                'DP_QUANT_VALUE2, DP_INTEGER2, DP_UNIT_CODE2,DP_CONFIDENTIAL2_YN, '\
                'CREATED_BY, CREATED_DATE from data_point')

        data_point = []
        for line in oracleCur:
            line = list(line)
            if line[9] != None:
                line[9] = line[9].read()
            data_point.append(line)

        self.localDb.executemany('insert into data_point values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',data_point)
        del data_point
            
        print 'data_point table done'

        oracleCur.execute('select QUALIFIER_ID, DATA_POINT_ID, Q_L0_DOMAIN_CODE, '\
                'Q_L1_KINGDOM_CODE, Q_QUAL_VALUE1, Q_QUANT_VALUE1, Q_INTEGER1, Q_UNIT_CODE1, '\
                'Q_CONFIDENTIAL1_YN, Q_CLOB1, Q_QUAL_VALUE2, Q_QUANT_VALUE2, Q_INTEGER2, '\
                'Q_UNIT_CODE2, Q_CONFIDENTIAL2_YN,CREATED_BY, CREATED_DATE from qualifier')

        qualifier = []
        for line in oracleCur:
            line = list(line)
            if line[9] != None:
                line[9] = line[9].read()
            qualifier.append(line)

        self.localDb.executemany('insert into qualifier values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',qualifier)
        del qualifier

        print 'qualifier done'

        self.localDb.execute('create index chemIdx on chemical (chemical_id)')
        print 'chemical index done'
        self.localDb.execute('create index dpIdx on data_point (chemical_id,data_point_id)')
        print 'data_point index done'
        self.localDb.execute('create index qIdx on qualifier (data_point_id,qualifier_id)')
        print 'qualifier index done'

        #sqlQuery =raw_input('Input sql query (die to quit)')
        #while sqlQuery.strip()!='die':
        #    try:
        #        results = list(localDb.execute(sqlQuery))
        #        for item in results:
        #            print ' '.join(item)
        #    except sqlite3.Error,dberror:
        #        print dberror
        #    sqlQuery = raw_input('SQLITE3>')
        #results = self.localDb.execute('select data_point_id from qualifier')
        #pprint(list(results))

    def cleanup(self):
        self.oracleCon.close()

    def __enter__(self):
        return self
    def __exit__(self,type,value,traceback):
        self.cleanup()
        
            

    def pullXml(self,saveFileName):
        #make root element
        
        with open(saveFileName,'w') as saveFile:
            saveFile.write('<ROOT>\n')
            #cur.execute('select * from data_point where chemical_id = 102')
            #pprint(cur.fetchall())
            #sys.exit()

            cur1 = self.localDb.cursor()
            cur2 = self.localDb.cursor()
            
            chemno = 0
            for chemid in self.chemids:
                chemno += 1
                print chemno
                #make parent chemical element
                chemical = etree.Element('CHEMICAL')
                dpchemical = self.localDb.execute('select chemical_id,created_by,created_date from chemical where chemical_id ='+str(chemid))

                chemelementlist = ['CHEMICAL_ID','CREATED_BY','CREATED_DATE']
                for c in dpchemical:
                    if len(c) != len(chemelementlist):
                        print 'len(c) != len(chemelementlist)'
                        sys.exit()
                        
                    for index in range(len(c)):
                        if str(c[index]) != 'None':
                            data = etree.Element(chemelementlist[index])
                            if chemelementlist[index] == 'CREATED_DATE' or chemelementlist[index] == 'MODIFIED_DATE':
                                #print c[index],type(c[index])
                                data.text = c[index].isoformat()
                                #Hack to fix
                            else:
                                data.text = str(c[index])
                            chemical.append(data)
                            
                
                #now get the datapoints.
                datapoints = etree.Element('datapoints')
    ##            dpResults = cur1.execute('select DATA_POINT_ID,CHEMICAL_ID,\
    ##                DP_L0_DOMAIN_CODE,DP_L1_KINGDOM_CODE,DP_QUANT_VALUE1,DP_QUANT_VALUE2,\
    ##                DP_QUAL_VALUE1,DP_QUAL_VALUE2,DP_CLOB1,DP_INTEGER1,DP_CONFIDENTIAL1_YN,\
    ##                DP_CONFIDENTIAL2_YN,CREATED_BY,CREATED_DATE,DP_UNIT_CODE1,DP_UNIT_CODE2 \
    ##                from data_point where chemical_id ='+str(chemid))

                dpResults = cur1.execute('select DATA_POINT_ID,CHEMICAL_ID,\
                    DP_L0_DOMAIN_CODE,DP_L1_KINGDOM_CODE,CREATED_BY,CREATED_DATE,\
                    DP_QUANT_VALUE1,DP_QUANT_VALUE2,\
                    DP_CLOB1,DP_INTEGER1,DP_CONFIDENTIAL1_YN,DP_CONFIDENTIAL2_YN,\
                    DP_QUAL_VALUE1,\
                    DP_QUAL_VALUE2,DP_UNIT_CODE1,DP_UNIT_CODE2\
                    from data_point where chemical_id ='+str(chemid))
                
                #dbdatapoints = self.cur.fetchall()
                #dbdatapoints.sort()
                #pprint(dbdatapoint)

                elementlist = ['DATA_POINT_ID','CHEMICAL_ID','DP_L0_DOMAIN_CODE','DP_L1_KINGDOM_CODE',\
                        'CREATED_BY','CREATED_DATE','DP_QUANT_VALUE1','DP_QUANT_VALUE2',\
                        'DP_CLOB1','DP_INTEGER1',\
                        'DP_CONFIDENTIAL1_YN','DP_CONFIDENTIAL2_YN',\
                        'DP_QUAL_VALUE1',\
                        'DP_QUAL_VALUE2','DP_UNIT_CODE1','DP_UNIT_CODE2']
                
                for dp in dpResults:
                    #print dp[2]
                    if len(dp) != len(elementlist):
                        print 'len(dp) != len(elementlist)'
                        sys.exit()
                    
                    #make datapoint tag
                    dpelement = etree.Element(dp[2])
    ##                for index in range(len(dp)):
    ##                    #rint elementlist[index]
    ##                    # omit DP_L0_DOMAIN_CODE - it's already the tag name of dpelement.
    ##                    if index != 2 and str(dp[index]) != 'None':
    ##                        data = etree.Element(elementlist[index])
    ##                        if elementlist[index] == 'CREATED_DATE':
    ##                            data.text = dp[index].isoformat()
    ##                        else:
    ##                            data.text = str(dp[index])
    ##                        if dp[2] == 'NICNASID' or dp[2] == 'CASNUMBER' and elementlist[index] == 'DP_QUANT_VALUE1':
    ##                            #fix integers - stored as floats in database. Remove .0
    ##                            #data.text = data.text.rstrip('0').rstrip('.')
    ##                            data.text = re.sub('\.0$','',data.text)
    ##                        dpelement.append(data)

                    for index in [0,1,3,4]:
                        if dp[index] != None:
                            data = etree.Element(elementlist[index])
                            data.text = str(dp[index])
                            dpelement.append(data)
                    if dp[5] != None:
                        data = etree.Element(elementlist[5])
                        
                        data.text=dp[5].isoformat()
                        dpelement.append(data)

                    for index in [6,7,8,9,10,11,12,13,14,15]:
                        if dp[index] != None:
                            data = etree.Element(elementlist[index])
                            data.text = str(dp[index])
                            dpelement.append(data)
                    
                    #now add the qualifiers.
                    dpid = dp[0]
                    #print dp[0],type(dp[0])
                    #qResults = cur2.execute('select QUALIFIER_ID,DATA_POINT_ID,Q_L0_DOMAIN_CODE,\
                    #    Q_L1_KINGDOM_CODE,Q_QUANT_VALUE1,Q_QUANT_VALUE2,Q_QUAL_VALUE1,Q_QUAL_VALUE2,\
                    #    Q_CLOB1,Q_INTEGER1,Q_CONFIDENTIAL1_YN,Q_CONFIDENTIAL2_YN,CREATED_BY,CREATED_DATE,\
                    #    Q_UNIT_CODE1,Q_UNIT_CODE2 from qualifier where data_point_id ='+str(dpid))
                    qResults = cur2.execute('select QUALIFIER_ID,DATA_POINT_ID,Q_L0_DOMAIN_CODE,\
                    Q_L1_KINGDOM_CODE,CREATED_BY,CREATED_DATE,Q_QUANT_VALUE1,Q_QUANT_VALUE2,\
                    Q_CLOB1,Q_INTEGER1,Q_CONFIDENTIAL1_YN,\
                    Q_CONFIDENTIAL2_YN,\
                    Q_QUAL_VALUE1,Q_QUAL_VALUE2,Q_UNIT_CODE1,Q_UNIT_CODE2 from qualifier where data_point_id ='+str(dpid))
                    #print list(qResults)

     

                    #dbqualifiers = self.cur.fetchall()
                    qualifiers = etree.Element('qualifiers')
                    qelementlist = ['QUALIFIER_ID','DATA_POINT_ID','Q_L0_DOMAIN_CODE',\
                            'Q_L1_KINGDOM_CODE','CREATED_BY',\
                            'CREATED_DATE','Q_QUANT_VALUE1','Q_QUANT_VALUE2',\
                            'Q_CLOB1','Q_INTEGER1',\
                            'Q_CONFIDENTIAL1_YN','Q_CONFIDENTIAL2_YN',\
                            'Q_QUAL_VALUE1','Q_QUAL_VALUE2',\
                            'Q_UNIT_CODE1','Q_UNIT_CODE2']
    ##                for q in qResults:
    ##                    if len(q) != len(qelementlist):
    ##                        print 'len(q) != len(qelementlist)'
    ##                        sys.exit()
    ##                    #make qualifier tag:
    ##                    qelement = etree.Element(q[2])
    ##                    for index in range(len(q)):
    ##                        if index != 2 and str(q[index]) != 'None':
    ##                            data = etree.Element(qelementlist[index])
    ##                            if qelementlist[index] == 'CREATED_DATE':
    ##                                data.text=q[index].isoformat()
    ##                            else:
    ##                                data.text = str(q[index])
    ##                            qelement.append(data)
    ##                    qualifiers.append(qelement)

                    for q in qResults:

                        #make qualifier tag (Unroll loops for speed).
                        qelement = etree.Element(q[2])
                        for index in [0,1,3,4]:
                            if q[index] != None:
                                data = etree.Element(qelementlist[index])
                                data.text = str(q[index])
                                qelement.append(data)
                        if q[5] != None:
                            data = etree.Element(qelementlist[5])
                            data.text=q[5].isoformat()
                            qelement.append(data)
                            
                        for index in [6,7,8,9,10,11,12,13,14,15]:
                            if q[index] != None:
                                data = etree.Element(qelementlist[index])
                                data.text = str(q[index])
                                qelement.append(data)
                        
                            
                        qualifiers.append(qelement)



                    for qualifier in qualifiers:
                        dpelement.append(qualifier)

                    #cur.execute('select * from qualifier where data_point_id ='+str(dpid))
                    #pprint(cur.fetchall())
                    #cur.execute('select column_name from all_tab_columns where table_name=\'QUALIFIER\'')
                    #pprint(cur.fetchall())
                    #sys.exit()
                    

                    datapoints.append(dpelement)
                    
                for datapoint in datapoints:
                    chemical.append(datapoint)

                #write chemical as soon as it is made
                try:
                    saveFile.write(etree.tostring(chemical,pretty_print=True))
                except IOError:
                    QtGui.QMessageBox.information(self,'Cascading Viewer',\
                        'Couldn\'t save for some reason.',\
                        QtGui.QMessageBox.Ok)
            saveFile.write('</ROOT>\n')
        print 'Xml download complete.'

        

