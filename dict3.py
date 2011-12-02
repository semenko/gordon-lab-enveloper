from lxml import etree


NAMES = "{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}"

infile = 'in.xml'
data = etree.parse(infile)

scans = data.findall("//%sscan" % NAMES)

for scan in scans:
   t print scan
