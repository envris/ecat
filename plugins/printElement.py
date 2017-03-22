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


def printFloat(fle):
    def isGt(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isLt(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' not in keys \
                and 'DP_QUANT_VALUE1' not in keys\
                and 'DP_QUAL_VALUE2' in keys\
                and 'DP_QUANT_VALUE2' in keys:
            return True
        else:
            return False
    
    def isRange(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' in keys\
                and 'DP_QUANT_VALUE2' in keys:
            return True
        else:
            return False

    def isSingle(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' not in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isExcept(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' not in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False




    #get data out of float element
    if fle == [] or fle == None:
        return '??'
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    keys = leaves.keys()
    #print leaves

    if 'DP_L1_KINGDOM_CODE' not in leaves:
        return '?'
    
    
    if isSingle(leaves) and 'DP_UNIT_CODE1' not in keys:
        return leaves['DP_QUANT_VALUE1']
    if isSingle(leaves) and 'DP_UNIT_CODE1' in keys:
        return leaves['DP_QUANT_VALUE1']+' '+leaves['DP_UNIT_CODE1']
    if isGt(leaves) and 'DP_UNIT_CODE1' not in keys:
        return leaves['DP_QUAL_VALUE1']+leaves['DP_QUANT_VALUE1']
    if isGt(leaves) and 'DP_UNIT_CODE1' in keys:
        return leaves['DP_QUAL_VALUE1']+leaves['DP_QUANT_VALUE1']+' '+leaves['DP_UNIT_CODE1']
    if isLt(leaves) and 'DP_UNIT_CODE2' not in keys:
        return leaves['DP_QUAL_VALUE2']+leaves['DP_QUANT_VALUE2']
    if isLt(leaves) and 'DP_UNIT_CODE2' in keys:
        return leaves['DP_QUAL_VALUE2']+leaves['DP_QUANT_VALUE2']+' '+leaves['DP_UNIT_CODE2']
    returnValue = ''
    if isRange(leaves):
        if leaves['DP_QUAL_VALUE1'] == '>':
            l = '('
        elif leaves['DP_QUAL_VALUE1'] == '>=':
            l = '['
        else:
            l=''
        
        if leaves['DP_QUAL_VALUE2'] == '<':
            r = ')'
        elif leaves['DP_QUAL_VALUE2'] == '<=':
            r = ']'
        else:
            r=''

        if 'DP_UNIT_CODE1' not in keys and  'DP_UNIT_CODE2' not in keys:
            return l+leaves['DP_QUANT_VALUE1']+','+leaves['DP_QUANT_VALUE2']+r
        if 'DP_UNIT_CODE1' in keys and  'DP_UNIT_CODE2' in keys:
            return l+leaves['DP_QUANT_VALUE1']+','+leaves['DP_QUANT_VALUE2']+r+' '+leaves['DP_UNIT_CODE1']
        return '?'
    if isExcept(leaves):
        return leaves['DP_QUAL_VALUE1']


def printVarchar2(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    return leaves['DP_QUAL_VALUE1']

def printYN(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    return leaves['DP_QUAL_VALUE1']


def printCLOB(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    return leaves['DP_CLOB1']


def printInteger(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    return leaves['DP_INTEGER1']


def printVarchar2Choice(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    return leaves['DP_QUAL_VALUE1']

def printDPElement(fle):
    def isGt(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isLt(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' not in keys \
                and 'DP_QUANT_VALUE1' not in keys\
                and 'DP_QUAL_VALUE2' in keys\
                and 'DP_QUANT_VALUE2' in keys:
            return True
        else:
            return False
    
    def isRange(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' in keys\
                and 'DP_QUANT_VALUE2' in keys:
            return True
        else:
            return False

    def isSingle(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' not in keys \
                and 'DP_QUANT_VALUE1' in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isExcept(leaves):
        keys = leaves.keys()
        if 'DP_QUAL_VALUE1' in keys \
                and 'DP_QUANT_VALUE1' not in keys\
                and 'DP_QUAL_VALUE2' not in keys\
                and 'DP_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}

    if leaves['DP_L1_KINGDOM_CODE'] == 'DP_YN' or \
            leaves['DP_L1_KINGDOM_CODE'] == 'DP_VARCHAR2'\
            or leaves['DP_L1_KINGDOM_CODE'] == 'DP_VARCHAR2CHOICE':
        return leaves['DP_QUAL_VALUE1']
    elif leaves['DP_L1_KINGDOM_CODE'] == 'DP_INTEGER':
        return leaves['DP_INTEGER1']
    elif leaves['DP_L1_KINGDOM_CODE'] == 'DP_CLOB':
        return leaves['DP_CLOB1']
    elif leaves['DP_L1_KINGDOM_CODE'] == 'DP_FLOAT':

        #get data out of float element
        if fle == [] or fle == None:
            return '??'
        leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
        keys = leaves.keys()
        #print leaves

        if 'DP_L1_KINGDOM_CODE' not in leaves:
            return '?'
        
        
        if isSingle(leaves) and 'DP_UNIT_CODE1' not in keys:
            return leaves['DP_QUANT_VALUE1']
        if isSingle(leaves) and 'DP_UNIT_CODE1' in keys:
            return leaves['DP_QUANT_VALUE1']+' '+leaves['DP_UNIT_CODE1']
        if isGt(leaves) and 'DP_UNIT_CODE1' not in keys:
            return leaves['DP_QUAL_VALUE1']+leaves['DP_QUANT_VALUE1']
        if isGt(leaves) and 'DP_UNIT_CODE1' in keys:
            return leaves['DP_QUAL_VALUE1']+leaves['DP_QUANT_VALUE1']+' '+leaves['DP_UNIT_CODE1']
        if isLt(leaves) and 'DP_UNIT_CODE2' not in keys:
            return leaves['DP_QUAL_VALUE2']+leaves['DP_QUANT_VALUE2']
        if isLt(leaves) and 'DP_UNIT_CODE2' in keys:
            return leaves['DP_QUAL_VALUE2']+leaves['DP_QUANT_VALUE2']+' '+leaves['DP_UNIT_CODE2']
        returnValue = ''
        if isRange(leaves):
            if leaves['DP_QUAL_VALUE1'] == '>':
                l = '('
            elif leaves['DP_QUAL_VALUE1'] == '>=':
                l = '['
            else:
                l=''
            
            if leaves['DP_QUAL_VALUE2'] == '<':
                r = ')'
            elif leaves['DP_QUAL_VALUE2'] == '<=':
                r = ']'
            else:
                r=''

            if 'DP_UNIT_CODE1' not in keys and  'DP_UNIT_CODE2' not in keys:
                return l+leaves['DP_QUANT_VALUE1']+','+leaves['DP_QUANT_VALUE2']+r
            if 'DP_UNIT_CODE1' in keys and  'DP_UNIT_CODE2' in keys:
                return l+leaves['DP_QUANT_VALUE1']+','+leaves['DP_QUANT_VALUE2']+r+' '+leaves['DP_UNIT_CODE1']
            return '?'
        if isExcept(leaves):
            return leaves['DP_QUAL_VALUE1']
    else:
        print 'printElement: Unknown element type '+leaves['DP_L1_KINGDOM_CODE']
        exit()

def printQElement(fle):
    def isGt(leaves):
        keys = leaves.keys()
        if 'Q_QUAL_VALUE1' in keys \
                and 'Q_QUANT_VALUE1' in keys\
                and 'Q_QUAL_VALUE2' not in keys\
                and 'Q_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isLt(leaves):
        keys = leaves.keys()
        if 'Q_QUAL_VALUE1' not in keys \
                and 'Q_QUANT_VALUE1' not in keys\
                and 'Q_QUAL_VALUE2' in keys\
                and 'Q_QUANT_VALUE2' in keys:
            return True
        else:
            return False
    
    def isRange(leaves):
        keys = leaves.keys()
        if 'Q_QUAL_VALUE1' in keys \
                and 'Q_QUANT_VALUE1' in keys\
                and 'Q_QUAL_VALUE2' in keys\
                and 'Q_QUANT_VALUE2' in keys:
            return True
        else:
            return False

    def isSingle(leaves):
        keys = leaves.keys()
        if 'Q_QUAL_VALUE1' not in keys \
                and 'Q_QUANT_VALUE1' in keys\
                and 'Q_QUAL_VALUE2' not in keys\
                and 'Q_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    def isExcept(leaves):
        keys = leaves.keys()
        if 'Q_QUAL_VALUE1' in keys \
                and 'Q_QUANT_VALUE1' not in keys\
                and 'Q_QUAL_VALUE2' not in keys\
                and 'Q_QUANT_VALUE2' not in keys:
            return True
        else:
            return False

    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}

    if leaves['Q_L1_KINGDOM_CODE'] == 'Q_YN' or \
            leaves['Q_L1_KINGDOM_CODE'] == 'Q_VARCHAR2'\
            or leaves['Q_L1_KINGDOM_CODE'] == 'Q_VARCHAR2CHOICE':
        return leaves['Q_QUAL_VALUE1']
    elif leaves['Q_L1_KINGDOM_CODE'] == 'Q_INTEGER':
        return leaves['Q_INTEGER1']
    elif leaves['Q_L1_KINGDOM_CODE'] == 'Q_CLOB':
        return leaves['Q_CLOB1']
    elif leaves['Q_L1_KINGDOM_CODE'] == 'Q_FLOAT':

        #get data out of float element
        if fle == [] or fle == None:
            return '??'
        leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
        keys = leaves.keys()
        #print leaves

        if 'Q_L1_KINGDOM_CODE' not in leaves:
            return '?'
        
        
        if isSingle(leaves) and 'Q_UNIT_CODE1' not in keys:
            return leaves['Q_QUANT_VALUE1']
        if isSingle(leaves) and 'Q_UNIT_CODE1' in keys:
            return leaves['Q_QUANT_VALUE1']+' '+leaves['Q_UNIT_CODE1']
        if isGt(leaves) and 'Q_UNIT_CODE1' not in keys:
            return leaves['Q_QUAL_VALUE1']+leaves['Q_QUANT_VALUE1']
        if isGt(leaves) and 'Q_UNIT_CODE1' in keys:
            return leaves['Q_QUAL_VALUE1']+leaves['Q_QUANT_VALUE1']+' '+leaves['Q_UNIT_CODE1']
        if isLt(leaves) and 'Q_UNIT_CODE2' not in keys:
            return leaves['Q_QUAL_VALUE2']+leaves['Q_QUANT_VALUE2']
        if isLt(leaves) and 'Q_UNIT_CODE2' in keys:
            return leaves['Q_QUAL_VALUE2']+leaves['Q_QUANT_VALUE2']+' '+leaves['Q_UNIT_CODE2']
        returnValue = ''
        if isRange(leaves):
            if leaves['Q_QUAL_VALUE1'] == '>':
                l = '('
            elif leaves['Q_QUAL_VALUE1'] == '>=':
                l = '['
            else:
                l=''
            
            if leaves['Q_QUAL_VALUE2'] == '<':
                r = ')'
            elif leaves['Q_QUAL_VALUE2'] == '<=':
                r = ']'
            else:
                r=''

            if 'Q_UNIT_CODE1' not in keys and  'Q_UNIT_CODE2' not in keys:
                return l+leaves['Q_QUANT_VALUE1']+','+leaves['Q_QUANT_VALUE2']+r
            if 'Q_UNIT_CODE1' in keys and  'Q_UNIT_CODE2' in keys:
                return l+leaves['Q_QUANT_VALUE1']+','+leaves['Q_QUANT_VALUE2']+r+' '+leaves['Q_UNIT_CODE1']
            return '?'
        if isExcept(leaves):
            return leaves['Q_QUAL_VALUE1']
    else:
        print 'printElement: Unknown element type '+leaves['Q_L1_KINGDOM_CODE']
        exit()

def printElement(fle):
    leaves= {str(x.tag): str(x.text) for x in fle.xpath('./*[not(*)]')}
    keys = leaves.keys()
    if 'DP_L1_KINGDOM_CODE' in keys:
        return printDPElement(fle)
    elif 'Q_L1_KINGDOM_CODE' in keys:
        return printQElement(fle)
    else:
        return None



