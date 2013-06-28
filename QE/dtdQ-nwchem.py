import re
import sys

CmtoEv = 0.0001239842573


output = open(sys.argv[1], 'r')

freqs_file = open(sys.argv[2], 'r')
freqs = freqs_file.readlines()
freqs_file.close()

norb = sys.argv[3]

i = 0
t = 0.0
L = 0.0
while True:
    line = output.readline()
    if not line:
        break
    p = re.search(r'^\s*'+norb+'\s*'+norb+'\s*(-*\d+\.\d+)', line)
    if p:
        if i % 2 == 0 and i > 0:
            if float(freqs[i/2-1]) == 0.0:
                print "Zero freq at " + str(i)
            else:
                v = abs(abs(t)-abs(float(p.group(1))))/0.01
                hw = float(freqs[i/2-1])*CmtoEv*1000.0 # meV
                L += v**2/(2*hw)
                print 'for w = %5.3f (cm-1): v = %5.3f ' % (float(freqs[i/2-1]), v)
        else:
            t = float(p.group(1))

        i += 1

output.close()
print 'L = %6.4f' % L