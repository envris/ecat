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


import sys
import os
from lxml import etree
import re
from copy import deepcopy


class makeTags(object):
    def __init__(self):
        super(makeTags,self).__init__()

    def makeFloat(self,dpTag,dpData,unit=None):
        #dpType = self.columnHeaderDict[colHeader]['type']
        #dpTag = self.columnHeaderDict[colHeader]['tag']
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
            
            #parse
            if re.match('^ *(?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?$',dpData):
                #it's a straight number.
    
                f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                f1.text = dpData
    
    
                confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                confidential.text = 'N'
    
                #unit = self.columnHeaderDict[colHeader]['unit']
                if unit != None:
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = unit
                #if re.match('^(?:<|<=|>|>=)? *(?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?$',dp):
                #pass
            elif re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData):
                #(>,>=,<,<=)
                gtlt = re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData).group(1)
                num = re.match('^(<|<=|>|>=)? *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData).group(2)
            
                #lower bound.
                if gtlt == '>' or gtlt == '>=':
                    f1 = etree.SubElement(Xml,'DP_QUANT_VALUE1')
                    f1.text = num
                    
    
                    confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
                    confidential.text = 'N'
        
                    q1 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
                    q1.text = gtlt
            
                    #units
                    #unit = self.columnHeaderDict[colHeader]['unit']
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
                    #unit = self.columnHeaderDict[colHeader]['unit']
                    if unit != None:
                        u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                        u2.text = unit
            
            elif re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?) *to *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData):
                # a to b
    
                atobmatch = re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?) *to *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData)                
    
    
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
                
                #unit = self.columnHeaderDict[colHeader]['unit']
                if unit != None:
                    u1 =  etree.SubElement(Xml,'DP_UNIT_CODE1')
                    u1.text = unit
                    u2 =  etree.SubElement(Xml,'DP_UNIT_CODE2')
                    u2.text = unit
    
            elif re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?) *- *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData):
                # a - b
    
                atobmatch = re.match('^ *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?) *- *((?:(?:-|\+)? *)(?:[0-9]+\.[0-9]+|\.[0-9]+|[0-9]+\.|[0-9]+)(?:(?:e|E)(?:\+|-)[0-9]+)?)$',dpData)
    
    
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
    
                #unit = self.columnHeaderDict[colHeader]['unit']
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
                errorReturn = [dpTag,dpData]
    
            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
    
    
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None
    
    def makeVarchar2(self,dpTag,dpData):
        #dpType = self.columnHeaderDict[colHeader]['type']
        #dpTag = self.columnHeaderDict[colHeader]['tag']
        #print colHeader, dpData,dpTag
        errorReturn = None
    
        if dpData != '' and dpData != '#N/A' and dpData != None:
    
            Xml = etree.Element(dpTag)
    
            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'
    
            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'
    
            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_VARCHAR2'
            
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'
    
            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
    
    
            confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
            confidential.text = 'N'
    
            varchar2 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
            varchar2.text = dpData
    
            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None

    
    def makeVarchar2Choice(self,dpTag,dpData):
        #dpType = self.columnHeaderDict[colHeader]['type']
        #dpTag = self.columnHeaderDict[colHeader]['tag']
        #print colHeader, dpData,dpTag
        errorReturn = None
    
        if dpData != '' and dpData != '#N/A' and dpData != None:
    
            Xml = etree.Element(dpTag)
    
            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'
    
            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'
    
            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_VARCHAR2CHOICE'
            
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'
    
            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
    
    
            confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
            confidential.text = 'N'
    
            varchar2 = etree.SubElement(Xml,'DP_QUAL_VALUE1')
            varchar2.text = dpData
    
            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None
    def makeClob(self,dpTag,dpData):
        errorReturn = None

        if dpData != '' and dpData != '#N/A' and dpData != None:

            Xml = etree.Element(dpTag)

            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'

            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'

            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_CLOB'

            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'A16922'

            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'

            clob = etree.SubElement(Xml,'DP_CLOB1')
            clob.text = dpData

            confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
            confidential.text = 'N'


            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None
  
    
    def makeYN(self,dpTag,dpData):
    
        #dpTag = self.columnHeaderDict[colHeader]['tag']
        #print colHeader, dpData,dpTag
        errorReturn = None
    
        if dpData != '' and dpData != '#N/A' and dpData != None:
    
            Xml = etree.Element(dpTag)
            
            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'
            
            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'
    
            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_YN'
            
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'
    
            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
    
            if dpData in ['Yes','1','Y','1.0']:
                dpData = 'Y'
            elif dpData in ['No','0','N','0.0']:
                dpData = 'N'
            elif dpData == 'Uncertain':
                dpData == 'Uncertain'
            elif dpData == '':
                return None
    
            else:
                #errorReturn = [self.data[index][self.fielddict['CAS-RN']],\
                #    colHeader,dpData]
                errorReturn = [colHeader,dpData]
            
    
            confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
            confidential.text = 'N'
    
            qv = etree.SubElement(Xml,'DP_QUAL_VALUE1')
            qv.text = dpData
    
    
            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None


    def makeInteger(self,dpTag,dpData):
        #dpType = self.columnHeaderDict[colHeader]['type']
        #dpTag = self.columnHeaderDict[colHeader]['tag']

        errorReturn = None
        
        if dpData != '' and dpData != '#N/A' and dpData != None:

            #parse.
            if re.match('^[0-9]+$',dpData):
                #it's an integer
                pass
            elif re.match('^[0-9]+\.0$',dpData):
                #it's an 'excel' integer. Remove trailing ".0"
                dpData = re.sub('\.0$','',dpData)
            else:
                #errorReturn = [self.data[index][self.fielddict['CAS-RN']],\
                #        colHeader,dpData]
                errorReturn = [colHeader,dpData]
            
            Xml = etree.Element(dpTag)

            dpid =  etree.SubElement(Xml,'DATA_POINT_ID')
            dpid.text = 'AUTO'

            chemid = etree.SubElement(Xml,'CHEMICAL_ID')
            chemid.text = 'AUTO'

            l1 = etree.SubElement(Xml,'DP_L1_KINGDOM_CODE')
            l1.text = 'DP_INTEGER'
            
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'

            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'
        
            i = etree.SubElement(Xml,'DP_INTEGER1')
            i.text = dpData
            
            confidential = etree.SubElement(Xml,'DP_CONFIDENTIAL1_YN')
            confidential.text = 'N'


            #for item in self.tagsDict[colHeader]:
            #        a = self.makeQTag(self.reportTagDict[item])
            #        Xml.append(a)
            if errorReturn != None:
                return errorReturn
            else:
                return Xml
        else:
            return None



    def makeQYN(self,name,YorN):
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

    def makeQTag(self,tagName):
        q = etree.Element(tagName)
        
        qid =  etree.SubElement(q,'QUALIFIER_ID')
        qid.text = 'AUTO'

        dpid =  etree.SubElement(q,'DATA_POINT_ID')
        dpid.text = 'AUTO'
        
        kc = etree.SubElement(q,'Q_L1_KINGDOM_CODE')
        kc.text = 'Q_TAG'
        
        cb = etree.SubElement(q,'CREATED_BY')
        cb.text = 'AUTO'
        
        cd = etree.SubElement(q,'CREATED_DATE')
        cd.text = 'AUTO'

        confidential = etree.SubElement(q,'Q_CONFIDENTIAL1_YN')
        confidential.text = 'N'

        return q

    def makeQClob(self,dpTag,dpData):
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
            l1.text = 'Q_CLOB'
            
            cb = etree.SubElement(Xml,'CREATED_BY')
            cb.text = 'AUTO'

            cd = etree.SubElement(Xml,'CREATED_DATE')
            cd.text = 'AUTO'

            clob = etree.SubElement(Xml,'Q_CLOB1')
            clob.text = dpData

            confidential = etree.SubElement(Xml,'Q_CONFIDENTIAL1_YN')
            confidential.text = 'N'


            return Xml
        else:
            return None

    def makeQFloat(self,qTag,data,unit):
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



    def makeQVarchar2(self,dpTag,dpData):
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
            
            clob = etree.SubElement(Xml,'Q_QUAL_VALUE1')
            clob.text = dpData
            


            return Xml
        else:
            return None

#if __name__ == '__main__':
#    print 'hello'
