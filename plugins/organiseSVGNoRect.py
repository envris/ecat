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
import glob
from pprint import pprint
from copy import deepcopy
import re

class svg:
    def __init__(self,svgString):
        
        self.currentx = None
        self.currenty = None

        self.bestx = None
        self.besty = None

        self.svg = etree.fromstring(svgString)


    def getWidth(self):
        internalSvg = self.svg.xpath('/svg:svg/svg:g/svg:svg',\
                namespaces={'svg': 'http://www.w3.org/2000/svg'})
        if len(internalSvg) == 1 and type(internalSvg) == list:
            viewBox = internalSvg[0].xpath('./@viewBox')
            if len(viewBox[0].split()) == 4:
                width = viewBox[0].split()[2]
                return float(width)
            else:
                return None
        else:
            return None

    def getHeight(self):
        internalSvg = self.svg.xpath('/svg:svg/svg:g/svg:svg',\
                namespaces={'svg': 'http://www.w3.org/2000/svg'})
        if len(internalSvg) == 1 and type(internalSvg) == list:
            viewBox = internalSvg[0].xpath('./@viewBox')
            if len(viewBox[0].split()) == 4:
                height = viewBox[0].split()[3]
                return float(height)
            else:
                return None
        else:
            return None

    def translateContents(self,x,y):
        g = etree.Element('g')
        if x == None:
            x = 0
        if y == None:
            y = 0
        g.set('transform', 'translate('+str(x)+','+str(y)+')')
        #for item in self.svg.xpath('/svg:svg/*',\
        #        namespaces={'svg': 'http://www.w3.org/2000/svg'}):
        #    g.append(item)

        ##get rid of title
        #for item in g.xpath('/svg:title',\
        #        namespaces={'svg': 'http://www.w3.org/2000/svg'}):
        #    item.getparent().remove(item)

        #debugging rectangle
        #drect = etree.Element('rect')
        #drect.set('x','0')
        #drect.set('y','0')
        #drect.set('width',str(self.getWidth()))
        #drect.set('height',str(self.getHeight()))
        #drect.set('stroke-width','2')
        #drect.set('fill','none')
        #g.append(drect)

        #get the 'internal' svg (with the viewbox)
        internalSvg = self.svg.xpath('/svg:svg/svg:g/svg:svg',\
                    namespaces={'svg': 'http://www.w3.org/2000/svg'})
        if len(internalSvg) == 1 :
    
            svgCopy = deepcopy(internalSvg[0])

            for item in svgCopy.xpath('./*'):
                #print item.tag
                if item.tag != '{http://www.w3.org/2000/svg}rect':
                    g.append(item)

            #print etree.tostring(g,pretty_print=True)
            #print '###################################'

            #print etree.tostring(g,pretty_print=True)
            return g
        else:
            return None


class rect:
    def __init__(self,totalWidth,totalHeight):
        
        self.widths = [totalWidth]
        self.heights = [totalHeight]
        self.occMat = [[False]]
        
        #tolerance to take care of floating point problems when
        #fitting svgs into vacant spots.
        self.tol = 0.5

    def findVacant(self,width,height):
        #print 'self.widths: ', self.widths
        #print 'self.heights: ', self.heights
        #pprint(self.occMat)


        
        for j in range(len(self.widths)):
            for i in range(len(self.heights)):
                if self.occMat[i][j] == False:
                    #Have a candidate. If we put the svg here,
                    #will it be inside the bounding rectangle
                    #(totalWidth,totalHeight) ? Find the boxes
                    #that will contain the svg.

                    w = 0
                    h = 0

                    widthOK = False
                    for j2 in range(j,len(self.widths)):
                        w += self.widths[j2]
                        
                        if width <= w + self.tol:
                            #success!
                            widthOK = True
                            break


                    heightOK = False
                    for i2 in range(i,len(self.heights)):
                        h += self.heights[i2]

                        if  height <= h + self.tol:
                            #success!
                            heightOK = True
                            break

                    if widthOK == True and heightOK == True:
                        #j,j2,i,i2 hold bounding indices. Check to make
                        #sure all cells are vacant
                        vacant = True
                        for j3 in range(j,j2+1):
                            for i3 in range(i,i2+1):
                                if self.occMat[i3][j3] == True:
                                    #occupied.
                                    vacant = False
                                if vacant == False:
                                    break
                            if vacant == False:
                                break

                        if vacant == True:
                            #We have found a vacant range of cells that 
                            #the svg will fit into.
                            return {'leftIdx' : j, 'rightIdx' : j2,\
                                    'topIdx' : i, 'bottomIdx' : i2}

        #If we get here without returning, we haven't found a vacant spot.
        return None
    def getBoundingRect(self):
        # get rightmost filled column
        for j in reversed(range(len(self.widths))):
            colEmpty = True
            for i in range(len(self.heights)):
                if self.occMat[i][j] == True:
                    colEmpty = False
                    break
            if colEmpty == False:
                break

        #calculate width
        boundWidth = sum(self.widths[:j+1])

        #get leftmost empty row
        for i in reversed(range(len(self.heights))):
            rowEmpty = True
            for j in range(len(self.widths)):
                if self.occMat[i][j] == False:
                    rowEmpty = False
                    break
            if rowEmpty == False:
                break
        boundHeight = sum(self.heights[:i+1])

        return (boundWidth,boundHeight)




        #calculate height.






    def split(self,vacantSpot,width,height):
        leftIdx = vacantSpot['leftIdx']
        rightIdx = vacantSpot['rightIdx']
        topIdx = vacantSpot['topIdx']
        bottomIdx = vacantSpot['bottomIdx']
        
        vacantSpotWidth = sum(self.widths[leftIdx:rightIdx+1])
        vacantSpotHeight = sum(self.heights[topIdx:bottomIdx+1])

        if width < vacantSpotWidth:
            #update widths
            self.widths[rightIdx] = self.widths[rightIdx] - (vacantSpotWidth - width)
            self.widths.insert(rightIdx+1,vacantSpotWidth - width)
            
            #update self.occMat - insert column
            for row in self.occMat:
                row.insert(rightIdx+1,row[rightIdx])
        #elif width == vacantSpotWidth:
        elif width >= vacantSpotWidth:
            #need this to take floating point errors into account
            pass
        #else:
        #    print 'width must be less than vacantSpotWidth.'
        #    exit()

        if height < vacantSpotHeight:
            #update heights
            self.heights[bottomIdx] = self.heights[bottomIdx] \
                    - (vacantSpotHeight - height)
            self.heights.insert(bottomIdx+1,vacantSpotHeight - height)

            #update self.occMat - insert row
            newRow = []
            for item in self.occMat[bottomIdx]:
                newRow.append(item)
            self.occMat.insert(bottomIdx+1,newRow)
        #elif height == vacantSpotHeight:
        elif height >= vacantSpotHeight:
            #need this to take floating point errors into account
            pass
        #else:
        #    print 'height must be less than vacantSpotHeight.'
        #    exit()



    def occupy(self,vacantSpot):

        leftIdx = vacantSpot['leftIdx']
        rightIdx = vacantSpot['rightIdx']
        topIdx = vacantSpot['topIdx']
        bottomIdx = vacantSpot['bottomIdx']

        for j in range(leftIdx,rightIdx+1):
            for i in range(topIdx,bottomIdx+1):
                self.occMat[i][j] = True


class packer:
    def __init__(self,svgStringList):
        
        self.svgList = []
        for item in svgStringList:
            svgItem = svg(item)
            self.svgList.append(svgItem)

        #sort by height, then width
        self.svgList.sort(key=lambda a: (-a.getHeight(),-a.getWidth()))
        #print 'heights: ', [a.getHeight() for a in self.svgList]
        #print 'widths: ', [a.getWidth() for a in self.svgList]
        #for item in self.svgList:
        #    print item.getHeight(),item.getWidth() 
        #exit()

        self.maxTotalWidth = 0
        self.maxTotalHeight = 0

        #max width and height of any single svg 
        #svgList[0] is tallest -sorted)
        self.maxWidth = 0
        self.maxHeight = self.svgList[0].getHeight()
        #print 'self.maxHeight: ',self.maxHeight
        #print'---'  
        for item in self.svgList:
            width = item.getWidth()
            height = item.getHeight()
            #print '(width,height): ',(width,height)

            self.maxTotalWidth += width
            self.maxTotalHeight += height
            if width > self.maxWidth:
                self.maxWidth = width
        #print'---'  

        self.minWidth = self.maxWidth
        self.minHeight =self.maxHeight
        for item in self.svgList:
            width = item.getWidth()
            height = item.getHeight()

            if width < self.minWidth:
                self.minWidth = width
            if height < self.minHeight:
                self.minHeight = height

        self.xStep = self.minWidth*0.33
        self.yStep = self.minHeight*0.33
        #self.xStep = 1
        #self.yStep = 1

        #add one to give a bit of extra width.
        #This is necessary in the case of floating point
        #errors in calculating the widths.
        self.currentWidth = self.maxTotalWidth + 1
        self.currentHeight = self.maxHeight + 1

        #make sure the first packing is accepted.
        self.bestWidth = self.currentWidth + 2
        self.bestHeight = self.currentHeight + 2


        counter = 0
        #print 'self.maxWidth',self.maxWidth
        #print 'self.maxHeight',self.maxHeight
        #print 'self.maxTotalWidth ',self.maxTotalWidth
        #print 'self.maxTotalHeight ',self.maxTotalHeight
        #print 'self.currentWidth ',self.currentWidth
        #print 'self.currentHeight ',self.currentHeight
        #print '---'
        while self.currentWidth >= self.maxWidth \
                and self.currentHeight <= self.maxTotalHeight:
            counter += 1
            #print '--------------------------------------'
            #print 'counter,self.currentWidth,self.currentHeight'
            #print counter,self.currentWidth,self.currentHeight
            #print 'self.bestWidth,self.bestHeight'
            #print self.bestWidth,self.bestHeight
            #for item in self.svgList:
            #    print item.bestx,item.besty
            allFit = True

            #for item in self.svgList:
            #    print etree.tostring(item.svg,pretty_print=True)
            
            r = rect(self.currentWidth,self.currentHeight)
            c2 = 0
            for svgItem in self.svgList:
                c2 +=1 
                w = svgItem.getWidth()
                h = svgItem.getHeight()
                #print '(w,h): ',w,h
                vacantSpot = r.findVacant(w,h)
                #print 'c2,vacantSpot: ',c2,vacantSpot
                if vacantSpot != None:
                    #print 'splitting vacant spot'
                    r.split(vacantSpot,w,h)
                    #print r.occMat
                    r.occupy(vacantSpot)
                    #print 'occupying ...'
                    #print r.occMat
                    #print 'widths: ',r.widths
                    #print 'heights',r.heights
                    svgItem.currentx = sum(r.widths[:vacantSpot['leftIdx']])
                    svgItem.currenty = sum(r.heights[:vacantSpot['topIdx']])
                else:
                    allFit = False
                    break
            #print 'allFit: ',allFit

            if allFit == True:
                #print 'allFit = True',self.xStep
                
                #set current width and height to width and height
                #of true bounding rectangle.
                (self.currentWidth,self.currentHeight) = r.getBoundingRect()

                #minimise area of circle that has all four corners
                #of the rectangle on its circumference
                #(Area = (pi/4)*(x**2 + y**2)
                #if self.currentWidth**2 + self.currentHeight**2 <\
                #        self.bestWidth**2 + self.bestHeight**2:
                #
                #minimise area of rectangle
                #if self.currentWidth*self.currentHeight < \
                #        self.bestWidth*self.bestHeight:

                #minimise weight function. Contains a contribution from
                #the area or the rectangle and the area of the circumcircle.
                a = 0.0
                if a*(self.currentWidth * self.currentHeight) \
                        +(1.0-a)*3.14/4.0*(self.currentWidth**2 + self.currentHeight**2) \
                        < a*(self.bestWidth * self.bestHeight) \
                        +(1.0-a)*3.14/4.0*(self.bestWidth**2 + self.bestHeight**2):
                    #outFile = open('temp_'+str(counter)+'.svg','w')
                    #outFile.write(etree.tostring(self.makeSVG('current'),\
                    #        pretty_print=True))
                    #outFile.close()
                #
                    #update everything.
                    self.bestWidth = self.currentWidth
                    self.bestHeight = self.currentHeight
                    #print counter,self.bestWidth,self.bestHeight,self.bestWidth*self.bestHeight
                    for svgItem in self.svgList:
                        svgItem.bestx = svgItem.currentx
                        svgItem.besty = svgItem.currenty
                    
                #make the bounding box thinner
                self.currentWidth -= self.xStep
            elif allFit == False:
                #print 'allFit = False',self.yStep
                #make the bounding box taller
                self.currentHeight += self.yStep
            else:
                print '???'
                exit()
        #print '############'
        #print 'self.currentWidth ',self.currentWidth
        #print 'self.currentHeight ',self.currentHeight
        #print 'self.bestWidth,self.bestHeight'
        #print self.bestWidth,self.bestHeight
        #print '#######'
    def makeSVG(self,bestOrCurrent,withRect):

        outputWidth = 200
        
        if bestOrCurrent == 'best':
            aspectRatio = self.bestWidth/self.bestHeight
        elif bestOrCurrent == 'current':
            aspectRatio = self.currentWidth/self.currentHeight
        else:
            print 'makeSVG: ???'
            exit()

        outputHeight = str(outputWidth/aspectRatio)
        outputWidth = str(outputWidth)

        svgString = '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:cml="http://www.xml-cml.org/schema"  x="0" y="0" font-family="sans-serif" stroke="rgb(0,0,0)" stroke-width="1"  stroke-linecap="round"></svg>'
        opSvg = etree.fromstring(svgString)
        opSvg.set('width',str(outputWidth)+'px')
        opSvg.set('height',str(outputHeight)+'px')
        #print 'opSvg: ',opSvg
        vbTxt = '0 0 '+' '+outputWidth+' '+outputHeight
        opSvg.set('viewBox',vbTxt)

        #background rectangle
        backRect = etree.Element('rect')
        backRect.set('x','0')
        backRect.set('y','0')
        backRect.set('fill','white')
        backRect.set('width',outputWidth)
        backRect.set('height',outputHeight)
        backRect.set('stroke-width','0')

        opSvg.append(backRect)
 

        g1 = etree.Element('g')
        g1.set('transform','translate(0,0)')

        opSvg.append(g1)

        innerSvg = etree.Element('svg')
        g1.append(innerSvg)

        innerSvg.set('width',outputWidth)
        innerSvg.set('height',outputHeight)
        innerSvg.set('x','0')
        innerSvg.set('y','0')
        innerSvg.set('font-family','sans-serif')
        innerSvg.set('stroke','rgb(0,0,0)')
        innerSvg.set('stroke-width','2')
        innerSvg.set('stroke-linecap','round')

        if withRect == True:
            drect = etree.Element('rect')
            drect.set('x','0')
            drect.set('y','0')
            drect.set('stroke-width','2')
            drect.set('fill','none')
            if bestOrCurrent == 'current':
                drect.set('width',str(self.currentWidth))
                drect.set('height',str(self.currentHeight))
            elif bestOrCurrent == 'best':
                drect.set('width',str(self.bestWidth))
                drect.set('height',str(self.bestHeight))
            else:
                print 'bestOrCurrent must be best or current.'
                exit()
            g2 = etree.Element('g')
            g2.append(drect)
            innerSvg.append(g2)



        for index,item in enumerate(self.svgList):
            #print etree.tostring(item.translateContents(item.bestx,item.besty),\
            #        pretty_print=True)
            if bestOrCurrent == 'best':
                outItem = item.translateContents(item.bestx,item.besty)
            elif bestOrCurrent == 'current':
                  outItem = item.translateContents(item.currentx,item.currenty)
            innerSvg.append(outItem)
        #print 'innerSvg: ',innerSvg

 
        if bestOrCurrent == 'best':
            #print 'best'
            vbTxt = '0 0 '+str(self.bestWidth)+' '+str(self.bestHeight)
            #print 'vbTxt ', vbTxt
            innerSvg.set('viewBox',vbTxt)
        elif bestOrCurrent == 'current':
            #print 'current'
            vbTxt = '0 0 '+str(self.currentWidth)+' '+str(self.currentHeight)
            innerSvg.set('viewBox',vbTxt)
        else:
            print 'bestOrCurrent must be either \'best\' or \'current\''
            exit()

        return opSvg
