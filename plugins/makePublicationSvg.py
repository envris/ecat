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



import organiseSVGNoRect
import re
import openbabel
from lxml import etree
from subprocess import Popen, PIPE
import base64
import sys
from PyQt4 import QtGui,QtCore,QtSvg

def smiles2svg(inSmiles,withRects):
    """Generate a svg structure from a SMILES string."""
    #print 'smiles2svg'

    #Deal with smiles made up of multiple distinct structures
    smilesList = inSmiles.split('.')

    svgList = []
    obconv = openbabel.OBConversion()
    obconv.SetInAndOutFormats('smi','svg')
    for smiles in smilesList:
        mol = openbabel.OBMol()
        obconv.ReadString(mol,smiles)
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
        p = organiseSVGNoRect.packer(svgList)
        #print 'packed'
        if withRects == True:
            outSvg = p.makeSVG('best',True)
        elif withRects == False:
            outSvg = p.makeSVG('best',False)
            #outSvg = p.makeSVG('best',True)
        else:
            print 'smiles2svg: ??'
            exit()


        outSvg =  etree.tostring(outSvg,pretty_print=True)

        #outSvg = '<?xml version="1.0"?>\n'+outSvg
        return outSvg
    else:
        print 'smiles2svg: ???'
        exit()

def svg2jpg(svgString):
    svgMatch = re.match('.*viewBox.*viewBox=\"[ \t]*([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]+([0-9]+\.?[0-9]*)[ \t]*\".*$',svgString,re.DOTALL)
    if svgMatch:
        svgWidth = float(svgMatch.group(3))
        svgHeight = float(svgMatch.group(4))
    else:
        print 'makeSvgs: no match for width/height.'
        exit()
    
    #svgWidth = 200
    #svgHeight = 166
    #svgWidth = 100
    #svgHeight = 83
    
    #rect = QtCore.QRect(0,0,svgWidth,svgHeight)
    
    
    svgSize = QtCore.QSize(svgWidth,svgHeight)
    myXmlStreamReader = QtCore.QXmlStreamReader(svgString)
    myRenderer = QtSvg.QSvgRenderer(myXmlStreamReader)
    
    myImage = QtGui.QImage(svgSize,QtGui.QImage.Format_ARGB32)
    myTempPainter = QtGui.QPainter(myImage)
    myTempPainter.fillRect(myImage.rect(),QtCore.Qt.white)
    myRenderer.render(myTempPainter)
    myTempPainter.end()
    
    buf = QtCore.QBuffer()
    buf.open(QtCore.QIODevice.ReadWrite)
    myImage.save(buf,"JPG")
    
    buf.seek(0)
    jpg = buf.readAll()
    return jpg




def makeHtml(altFileName,svg):
    #svg = smiles2svg(smiles)

    svgb64 = base64.b64encode(svg)
    jpgb64 = base64.b64encode(svg2jpg(svg))

    html = '<object width="500" height="100%" type="image/svg+xml" data="data:image/svg+xml;base64,' +svgb64.strip() +  '"><img src="data:image/jpg;base64,'+ jpgb64.strip()\
            +'" alt="'+altFileName+'"/></object>'
    return html



def combinedSvgFromSmiles(smilesList):
    svgList = []
    for smiles in smilesList:
        svg = smiles2svg(smiles,False)
        svgList.append(svg)
    p = organiseSVGNoRect.packer(svgList)
    outSvg = p.makeSVG('best',False)

    #remove width and height attributes from top-level svg.
    del outSvg.attrib['width']
    del outSvg.attrib['height']


    outSvg= etree.tostring(outSvg,pretty_print=True)
    #print '#####################'
    #print outSvg

    #outString = makeHtmlFromSvg(re.sub('-','',key)+'jpg',outSvg)

    return outSvg


#if __name__=='__main__':
#    app = QtGui.QApplication(sys.argv,False)
#    
#    svg = combinedSvgFromSmiles(['[Ca++].CCCCCCCC(=O)[O-].CCCCCCCC(=O)[O-]'])
#    
#    print makeHtml('foo',svg,svg2jpg(svg))
