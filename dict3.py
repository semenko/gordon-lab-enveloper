from lxml import etree


NAMES = "{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}"

infile = 'in.xml'
data = etree.parse(infile)

scans = data.findall("//%sscan" % NAMES)

ms1 = {}
ms2 = {}

# {'polarity': '+', 'basePeakIntensity': '1.6726954e07', 'scanType': 'FULL', 'retentionTime': 'PT928.853S', 'basePeakMz': '875.460632324219', 'peaksCount': '14488', 'msLevel': '1', 'lowMz': '400.001117585287', 'num': '2825', 'scanEvent': '1', 'totIonCurrent': '1.31472064e08', 'highMz': '1814.504541581188', 'centroided': '0', 'msInstrumentID': 'IC1'}

ms1keys = ['polarity', 'basePeakIntensity', 'scanType', 'retentionTime', 'basePeakMz', 'peaksCount', 'lowMz', 'scanEvent', 'totIonCurrent', 'highMz', 'centroided']
ms1types = [str, float, str, str, float, int, float, int, float, float, int]


# {'polarity': '+', 'basePeakIntensity': '67.274475097656', 'scanType': 'FULL', 'collisionEnergy': '35.0', 'retentionTime': 'PT266.898S', 'basePeakMz': '633.969665527344', 'peaksCount': '26', 'msLevel': '2', 'lowMz': '363.079956054688', 'num': '810', 'scanEvent': '2', 'totIonCurrent': '235.690582275391', 'highMz': '1169.914184570313', 'centroided': '1', 'msInstrumentID': 'IC2'}
ms2keys = ['polarity', 'basePeakIntensity', 'scanType', 'collisionEnergy', 'retentionTime', 'basePeakMz', 'peaksCount', 'lowMz', 'scanEvent', 'totIonCurrent', 'highMz', 'centroided']
ms2types = [str, float, str, float, str, float, int, float, int, float, float, int]

# {'precursorIntensity': '1.9346121875e05', 'activationMethod': 'CID', 'precursorCharge': '1', 'precursorScanNum': '6993'}
# {'compressedLen': '0', 'pairOrder': 'm/z-int', 'precision': '32', 'byteOrder': 'network'}


for scan in scans:
#   print scan
#   print scan.tag
#   print scan.text
   print scan.attrib
   if scan.attrib['msInstrumentID'] == "IC1":
      # MS1
      ms1_temp_dict = {}
      for key, cast in zip(ms1keys, ms1types):
         ms1_temp_dict[key] = cast(scan.attrib[key])

      ms1[scan.attrib['num']] = ms1_temp_dict

      # 
      for child in scan:
         print child.items()
         print child.keys()
         for child2 in child:
            print child2

   else:
      # MS2
      ms2_temp_dict = {}
      for key, cast in zip(ms2keys, ms2types):
         ms2_temp_dict[key] = cast(scan.attrib[key])

      ms2[scan.attrib['num']] = ms2_temp_dict

#   print scan.attrib
#   print scan.attrib['num']
#   print scan.items()
#   print scan.keys()


   
   # These are the MS2 scans
   for child in scan:
      print child.attrib
   print "****"
