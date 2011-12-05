from lxml import etree


NAMES = "{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}"

infile = 'in.xml'
data = etree.parse(infile)

scans = data.findall("//%sscan" % NAMES)

for scan in scans:
   print scan
#   print scan.tag
#   print scan.text
   print scan.attrib
#   print scan.items()
#   print scan.keys()
   for child in scan:
      print child
   print "****"
