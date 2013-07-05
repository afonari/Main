import re
import sys
from math import sqrt
import numpy as np


def sym_mat(mat, size, debug=False):
    m = mat
    for i in range(size):
        for j in range(size):
            #dif = dif + abs(dynout(nu_i,nu_j)-dynout(nu_j,nu_i))
            m[i][j] = 0.5*(m[i][j]+m[j][i])
            m[j][i] = m[i][j]

    if debug == True:
        print "dif"

    return m

CmtoEv = 0.0001239842573
Bohr = 0.529177208590000 # [A]

scf_in = open(sys.argv[1], 'r')

stepsize = 0.005

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
#print symbols, masses, positions

dynmat_out = open(sys.argv[2], 'r')
q_point = []
eigvals = []
norms = []
eigenvecs = []

while True:
    line = dynmat_out.readline()
    if not line:
        break

    p = re.search(r'\s*q\s*=\s*([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)', line)
    if p:
        q_point = [float(p.group(1)), float(p.group(2)), float(p.group(3))]
        if any(x != 0.0 for x in q_point):
            print "Found q point which is not G point."
            sys.exit(0)

    p = re.search(r'\s*omega.+?([-\d\.]+)\s*\[cm-1\]; norm=\s*([-\d\.]+)', line)
    if p:
        eigval = float(p.group(1))
        norm = float(p.group(2))

        if eigval < 0.0:
            print "Skipping negative eigenvalue %5.3f cm-1" % eigval
            continue

        eigvals.append(eigval)
        norms.append(norm)
        for i in range(natoms):
            line = dynmat_out.readline()
            p1 = re.search(r'(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', line)
            if not p1:
                break

            eigenvec = [ (float(p1.group(1))+1j*float(p1.group(2))).real, (float(p1.group(3))+1j*float(p1.group(4))).real, (float(p1.group(5))+1j*float(p1.group(6))).real ]
            eigenvecs.append( eigenvec )

#stepsize = 0.005
coeff = 0.0
amu_ry = 911.444242132273
convfact = Bohr**2*sqrt(amu_ry)
omega_o_fpi = 116.023953513992

for i in range(len(eigvals)):

    ra_np = np.zeros((3,3))
    for step in (-1,1):
        if step == -1:
            coeff = -0.5
        else:
            coeff = 0.5

        filename = '%03d-%d' % (i, step)

        ph_file = open(filename+'-ph.in.out', 'r')

        while True:
            line = ph_file.readline()
            if not line:
                break

            p = re.search(r'Dielectric constant in cartesian axis', line)
            if p:
                ph_file.readline() # empty line
                ra =[]
                for j in range(3):
                    line = ph_file.readline()
                    p1 = re.search(r'\(\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', line)
                    line = [ float(p1.group(1)), float(p1.group(2)), float(p1.group(3)) ]
                    ra.append( line )

                #print ra
                ra_sym = sym_mat(ra, 3)
                #print ra_sym
                ra_np += np.array(ra_sym) * coeff/stepsize * sqrt(norms[i]) * omega_o_fpi * convfact

                break # while true

        if step == 1:
            alpha = (ra_np[0,0] + ra_np[1,1] + ra_np[2,2])/3.0
            beta2 = ( (ra_np[0,0] - ra_np[1,1])**2 + (ra_np[0,0] - ra_np[2,2])**2 + (ra_np[1,1] - ra_np[2,2])**2 + 6.0 * (ra_np[0,1]**2 + ra_np[0,2]**2 + ra_np[1,2]**2) )/2.0
            print eigvals[i], "  ", (45.0*alpha**2 + 7.0*beta2)#, alpha, beta2
            #print '%7.4f %7.4f %7.4f %7.4f %7.4f' % (ra_np[0,0], ra_np[0,1], ra_np[0,2], ra_np[1,1], ra_np[1,2])
            #print ra_np








