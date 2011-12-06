#from lxml import etree
from xml.etree import ElementTree as etree

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
precursorKeys = ['precursorIntensity', 'activationMethod', 'precursorScanNum']
precursorTypes = [float, str, int] # + data

peakKeys = ['compressedLen', 'pairOrder', 'precision', 'byteOrder']
peakTypes = [int, str, int, str] # + data

for scan in scans:
#   print scan
#   print scan.tag
#   print scan.text
#   print scan.attrib
   if scan.attrib['msInstrumentID'] == "IC1":
      # Store MS1 values in a dict
      ms1_temp_dict = {}
      for key, cast in zip(ms1keys, ms1types):
         ms1_temp_dict[key] = cast(scan.attrib[key])

      peak_temp_dict = {}
      # MS1 scans have only one peak child ("scan" is an iterable)
      # TODO: Remove these? Are these values useful?
      for key, cast in zip(peakKeys, peakTypes):
         peak_temp_dict[key] = cast(scan[0].attrib[key])
      peak_temp_dict['rawPeak'] = scan[0].text # The raw b64 data of the peak.
      ms1_temp_dict['peak'] = peak_temp_dict

      # Add them all to the final MS1 dict.
      ms1[scan.attrib['num']] = ms1_temp_dict

   else:
      # Store MS2 values in a dict
      if int(scan.attrib['peaksCount']) == 0:
         print "No peaks!"
         continue # Nothing to see here. Next iteration please.
      ms2_temp_dict = {}
      for key, cast in zip(ms2keys, ms2types):
         ms2_temp_dict[key] = cast(scan.attrib[key])

      peak_temp_dict = {}
      # MS2 scans have both "precursor" and "peak" children.
#      print scan[0].attrib
      for key, cast in zip(precursorKeys, precursorTypes):
         peak_temp_dict[key] = cast(scan[0].attrib[key])
      peak_temp_dict['precursorMz'] = scan[0].text # The raw precursor Mz

      for key, cast in zip(peakKeys, peakTypes):
         peak_temp_dict[key] = cast(scan[1].attrib[key])
      peak_temp_dict['rawPeak'] = scan[1].text # The raw b64 of the peak

      ms2_temp_dict['peak'] = peak_temp_dict

      ms2[scan.attrib['num']] = ms2_temp_dict


print len(ms1)
print len(ms2)
