#!/usr/bin/env python

from time import sleep

import xml.etree.ElementTree as etree

# From http://stackoverflow.com/questions/127606/editing-xml-as-a-dictionary-in-python

tree = etree.parse('in.xml')
root = tree.getroot()

def xml_to_dict(el):
  d={}
  if el.text:
    d[el.tag] = el.text
  else:
    d[el.tag] = {}
  children = el.getchildren()
  if children:
    d[el.tag] = map(xml_to_dict, children)
  return d

a = xml_to_dict(root)

mzXML = a['{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}mzXML'][0]['{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}msRun']

for element in mzXML:
  # Just in case some weird nesting takes place in a later schema or odd MS run type.
  assert(len(element) == 1)

  # Pick out just the scans (skip machine info, software version, etc.)
  if '{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}scan' in element:
    print element

sleep(1)
