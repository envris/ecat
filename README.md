ECAT - the Electronic Chemical Assessment Tool
==============================================

Introduction
------------
The Electronic Chemical Assessment Tool (ECAT) is a database application
designed by the Australian federal Department of the Environment and Energy
to aid the environmental assessment of chemicals under its existing chemicals 
programme.  Assessors use ECAT to enter chemical data, including structural 
information, ecotoxicological endpoints and physicochemical properties, store 
this data in a database and generate a variety of reports.  ECAT provides a 
powerful and intuitive searching GUI, allowing assessors to find and group 
chemicals based on their structure and properties. 

ECAT is written in python with a PyQt GUI and uses openbabel to perform
chemoinformatics operations.

ECAT is released under the GPL v2 licence.


Installation
------------
ECAT requires python 2.7 and the following python packages:

* PyQt4 (https://sourceforge.net/projects/pyqt/files/PyQt4/PyQt-4.11.3/. 
This version of PyQt is released under the GPL v2 licence.)
* openbabel (http://openbabel.org/docs/current/UseTheLibrary/PythonInstall.html)
* lxml (http://lxml.de/)

Some of the report plugins require bootstrap and jquery. Download the compressed
bootstrap javascript from getbootstrap.com and save it in the plugins 
directory as bootstrap.min.js. Download the cerulean theme for bootstrap 
from bootswatch.com and save in the plugins directory as bootstrap.min.css. 
Download the compressed jquery javascript from jquery.com and save in the 
plugins directory as jquery.min.js.

To use ECAT with a database backend, a database interface will also be required.
An example interface using cx_Oracle is found in dbXmlExample.py.

For further information, please contact CASadmin@environment.gov.au.

