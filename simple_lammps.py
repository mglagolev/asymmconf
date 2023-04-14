#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 16:01:52 2023

@author: Mikhail Glagolev
"""
def atom_string(atomid, molid, atomtype, x, y, z, ix, iy, iz):
	return " ".join([str(atomid).rjust(7),
                  str(molid).rjust(7),
                  str(atomtype).rjust(7),
                  ("%.03f" % x).rjust(17),
                  ("%.03f" % y).rjust(17),
                  ("%.03f" % z).rjust(17),
                  #str(x).rjust(17), str(y).rjust(17), str(z).rjust(17), 
                  str(ix).rjust(7), str(iy).rjust(7), str(iz).rjust(7)])

def vel_string(atomid, vx, vy, vz):
	return " ".join([str(atomid).rjust(7),
                  str(vx).rjust(17), str(vy).rjust(17), str(vz).rjust(17)])

def bond_string(bondid, bondtype, iatom1, iatom2):
	return " ".join([str(bondid).rjust(7), 
                  str(bondtype).rjust(7),
                  str(iatom1).rjust(7), str(iatom2).rjust(7)])

def angle_string(angleid, angletype, iatom1, iatom2, iatom3):
	return " ".join([str(angleid).rjust(7),
                  str(angletype).rjust(7), str(iatom1).rjust(7),
                  str(iatom2).rjust(7), str(iatom3).rjust(7)])

def dihedral_string(dihedralid, dihedraltype, iatom1, iatom2, iatom3, iatom4):
	return " ".join([str(dihedralid).rjust(7),
                  str(dihedraltype).rjust(7),
                  str(iatom1).rjust(7), str(iatom2).rjust(7),
                  str(iatom3).rjust(7), str(iatom4).rjust(7)])
