
import struct
import timeit
import base64


mystr = "Q0IK70BLQhNDRwxUQmvcyUNXH55AMK1PQ2dNgkCyGf9DezEMQQGAB0N8PORAxxVNQ4EGyELgyCVDhoNVQxCQv0OHg5RC5XanQ4vzBkAxEzBDj59+QCNkAEOQjbNDGn2TQ5Ku60CIsLRDk36hQRa9+EOaCGBC5xIAQ5qO9D+7eaBDmwjyQCNY10ObnFBBglS0Q5wgmkIXEzZDnIc6QpMWEUOjhO5AIzyuQ6QFdEAHH5RDpQIERTcFkUOlmRJBd4GsQ6aQxEE2LOBDqJPCQGjbFEOsFwRAmApYQ62KckAjZFhDsb02QGijdkPFjQ5BTm5UQ8hzikDooWRDynaYQIJkf0PPBFZAv9kmQ9B+ekBM7BFD0glQQJdhHkPS+UxAIvy3Q9kJiECrhuJD2vsmQD3U6EPcCuBATVLdQ98V2kDQJE1D4KZkQCn3RkPhJnJAjzkTQ+Z+0EAjOsxD5vs2QCK7O0Pnm4JBBn79Q+iT7kDMaYRD6YHoQCNVqEPqGG5BoFoGQ+qikECeHa9D7dXuQCNyK0PuS6hAMqtGQ/Cb0EAwrUJD9axgQAeDHkP8j3BA/lP1Q/4MikCwtaVD/p+oQAfET0QC8dhBbwVeRAMzCUDnsipEBFA+QBaVx0QExn9AprKARAXE7kNeeeZEBkLvQlA5HEQI3BtAWtZdRAlQOUCQNvhECcZsQomAvEQKA4ZBE1yuRApN6kDG9EBEC8RDQLLlGUQMB/xBJxTNRA2dUUCWcSREDjwOQMwigkQOs+E/8vd4RBBdk0Ifr3REEJk3QkvD3kQQxhBDNTb8RBEESUAjJLJEFL5HQIh0eA=="

decoded = base64.b64decode(mystr)
#decoded2 = base64.standard_b64decode(mystr)


#a = "ABC"
#print len(a)
#b = struct.pack('<3s', a)
#print len(b)
#print struct.unpack('>3s', b)

#print len(mystr)
#print len(decoded)
#print len(decoded)

#for i in range(0, len(decoded)/4, 2):
#    print i
#    print decoded[i*4:i*4+4]
#    print struct.unpack('!f', decoded[i*4:i*4+4])

#print len(range(0, (len(decoded)/4-1)*4, 8))

#for n in range(500000):
#    a = [struct.unpack('!ff', decoded[x:x+8]) for x in range(0, (len(decoded)/4-1)*4, 8)]

a = struct.unpack('>154f', decoded)
#print [(a[x], a[x+1]) for x in range(0, len(a), 2)]
print zip(a[::2], a[1::2])

#for n in range(10000000):
#    a = '>%if' % (len(decoded)/4)
#print '>' + str(len(decoded)/4) + 'f'

#for n in range(100000):
#    a = struct.unpack('>' + str(len(decoded)/4) + 'f', decoded)
 #   magic = zip(a[::2], a[1::2])
   

