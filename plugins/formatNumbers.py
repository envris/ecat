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

import re

#dp = decimal places
#function for HTML scientific notation
def HTMLscientific(num,dp):
	num = float(num)
	strnum = ('%.'+str(dp)+'e') % num

	[coefficient,exponent] = strnum.split('e')
	coefficient.strip()
	if re.search('-',exponent):
		exponent = exponent.lstrip('-')
		exponent = exponent.lstrip('0')
		exponent = '-'+exponent
	else:
		exponent = exponent.lstrip('+')
		exponent = exponent.lstrip('0')

	return coefficient + ' x 10<sup>' + exponent + '</sup>'

#function for LaTex scientific notation
def texscinot(num,dp):
	num = float(num)
	strnum = ('%.'+str(dp)+'e') % num

	[coefficient,exponent] = strnum.split('e')
	coefficient.strip()
	if re.search('-',exponent):
		exponent = exponent.lstrip('-')
		exponent = exponent.lstrip('0')
		exponent = '-'+exponent
	else:
		exponent = exponent.lstrip('+')
		exponent = exponent.lstrip('0')

	return '$'+coefficient + '\\times 10^{'+exponent+'}$'

#returns the number of significant figures
def getnumsf(x):
	if type(x)!=str:
		print 'Input to getnumsf needs to be a string.'
		exit()
	else:
		x=x.strip('-')
		x=x.strip()
		if (float(x)>=1) and not re.search('\.',x):
			sfstring = x.strip('0')
			return len(sfstring)
		elif (float(x)>=1) and re.search('\.',x):
			x.lstrip('0')
			sfstring = x.replace('.','')
			return len(sfstring)
		elif (float(x)<1):
			sfstring = x.lstrip('0')
			sfstring = sfstring.lstrip('\.')
			sfstring = sfstring.lstrip('0')
			return len(sfstring)
		else:
			print 'getnumsf: ???'
			exit()

def limitsf(limsf,item):

	if item.rstrip('0').rstrip('.').rstrip('0') == '': # if item is zero
		return '0'
	else:
		sf = getnumsf(item)
		outsf = sf
		if sf>limsf:
			outsf = limsf
		return ('%.'+str(outsf - 1)+'e') % float(item)

def printnum(item):
	try:
		nsf = getnumsf(item)
		if nsf == 0: # it's zero!
			return '0'
		else:
			num = float(limitsf(3,item))
			absnum = abs(num)
		
			if nsf >=4:
				nsf = 3
			
			#strnum = ('%.'+str(dp)+'e') % num
			if  absnum >= 0.0001 and absnum <0.001:
				dp = nsf + 3
			if  absnum >= 0.001 and absnum <0.01:
				dp = nsf + 2
			elif absnum >= 0.01 and absnum <0.1:
				dp = nsf + 1
			elif absnum >= 0.1 and absnum < 1.0:
				dp = nsf
			elif absnum >= 1.0 and absnum < 10.0:
				dp = nsf - 1
			elif absnum >= 10.0 and absnum < 100.0:
				dp = nsf -2
			elif absnum >= 100.0 and absnum < 1000.0:
				dp = nsf -3
			elif absnum >= 1000.0 and absnum < 10000.0:
				dp = nsf -4
			elif absnum >= 10000.0 and absnum < 100000.0:
				dp = nsf -5
			
			if absnum >= 0.0001 and absnum < 100000.0:
				if dp < 0:
					dp = 0
				if dp <= 5:
					strnum = ('%.'+str(dp)+'f') % num
				else:
					#strnum = ('%.'+str(nsf - 1)+'e') % num
					strnum = HTMLscientific(str(num),nsf-1)
			else:
				strnum = HTMLscientific(str(num),nsf-1)
			#return item+', '+strnum
			return strnum

	except:
		return item

def print2dp(value):
	try:
		 value = "%.2f" % float(value)
	except:
	         pass
	return value



