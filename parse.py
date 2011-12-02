#!/usr/bin/env python


import xml.etree.ElementTree as ET

tree = ET.parse('in.xml')

print tree

root = tree.getroot() #.findall('.//{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Internal-Data')[0].text

print "root"
print root

print root.attrib

element = root.findall('.//{http://sashimi.sourceforge.net/schema_revision/mzXML_3.1}mzXML')

print "element"
print element

print element.attrib
print element.tag
print element.text
print element.tail
print element.items()
print element.keys()
