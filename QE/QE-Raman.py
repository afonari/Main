import re
import sys

CmtoEv = 0.0001239842573

scf_in = open(sys.argv[1], 'r')

masses_dict = {}
symbols = []
masses = []
positions = []
while True:
    line = scf_in.readline()
    if not line:
        break

    p = re.search(r'ATOMIC_SPECIES', line)
    if p:
        while True:
            line = scf_in.readline()
            p1 = re.search(r'^\s*(\w+)\s*(\d+\.\d+)', line)

            if not p1:
                break
            masses_dict[p1.group(1)] = float(p1.group(2))

    p = re.search(r'ATOMIC_POSITIONS angstrom', line)
    if p:
        while True:
            line = scf_in.readline()
            p1 = re.search(r'^\s*(\w+)\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', line)

            if not p1:
                break

            atomic_symbol = p1.group(1)
            symbols.append( atomic_symbol )
            masses.append( masses_dict[atomic_symbol] )
            positions.append([ float(p1.group(2)), float(p1.group(3)), float(p1.group(4)) ])
scf_in.close()
natoms = len(symbols)
# print symbols, masses, positions

dynmat_out = open(sys.argv[2], 'r')
q_point = []
eigvals = []
eigenvecs = []

while True:
    line = dynmat_out.readline()
    if not line:
        break

    p = re.search(r'\s*q\s*=\s*\(?\s*([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)', line)
    if p:
        q_point = [float(p.group(1)), float(p.group(2)), float(p.group(3))]
        if any(x != 0.0 for x in q_point):
            print "Found q point which is not G point."
            sys.exit(0)

    p = re.search(r'\s*omega.+?([-\d\.]+)\s*\[cm-1\]', line)
    if p:
        eigval = float(p.group(1))
        if eigval < 0.0:
            print "Skipping negative eigenvalue %5.3f cm-1" % eigval
            continue

        eigvals.append(eigval)
        for i in range(natoms):
            line = dynmat_out.readline()
            p1 = re.search(r'(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', line)
            if not p1:
                break

            eigenvec = [ (float(p1.group(1))+1j*float(p1.group(2))).real, (float(p1.group(3))+1j*float(p1.group(4))).real, (float(p1.group(5))+1j*float(p1.group(6))).real ]
            eigenvecs.append( eigenvec )

stepsize = 0.01
for i in range(len(eigvals)):
    for step in (-1,1):
        filename = '%03d-%d' % (i, step)

        scf_file = open(filename+'-scf.in', 'w')
        ph_file = open(filename+'-ph.in', 'w')

        scf_in = open(sys.argv[1], 'r')
        while True:
            line = scf_in.readline()
            if not line:
                break

            p = re.search(r'\s*&control', line)
            if p:
                scf_file.write(line+'\n'+'    outdir = "./'+filename+'"\n')
                continue
            
            p = re.search(r'ATOMIC_POSITIONS angstrom', line)
            if p:
                for i in range(natoms):
                    pos = positions[i]
                    eigenvec = eigenvecs[i]
                    pos = [a + b for (a,b) in zip(pos,eigenvec)]
                    scf_file.write("%s %7.5f %7.5f %7.5f\n" % (symbols[i], pos[0], pos[1], pos[2]) )
                continue

            scf_file.write(line)






#print eigvals
#print eigenvecs

