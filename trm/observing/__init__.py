#!/usr/bin/env python

"""
Module to support observing-related routines such as eclipse predicting.
Provides some classes for reading in ephemeris and positional files and
a routine specific to the Thai telescope
"""

import math as m
import numpy as np
from trm import subs

class Ephemeris (object):
    """
    Stores an ephemeris.
    """
    def __init__(self, string):
        """
        Constructs an Ephemeris when passed a string having the form:

          Time Poly T0 eT0 P eP [Q eQ]

        Time    -- a time type: 'BMJD', 'HJD', or 'HMJD'
        Poly    -- 'linear', 'quadratic'
        T0, eT0 -- constant term and uncertainty
        P, eP   -- linear coefficient and uncertainty
        Q, eQ   -- quadtraic term and uncertainty
        """

        subv = string.split()

        self.time   = subv[0]
        if self.time != 'HJD' and self.time != 'BJD' and \
           self.time != 'HMJD' and self.time != 'BMJD':
            raise Exception('Ephem: unrecognised time type = ' + self.time)

        self.poly   = subv[1]
        if self.poly == 'linear':
            if len(subv) != 6:
                raise Exception('Ephem: linear ephemerides require 4 numbers')
        elif self.poly == 'quadratic':
            if len(subv) != 8:
                raise Exception('Ephem: quadratic ephemerides require 6 numbers')
        else:
            raise Exception("Ephem: only 'linear' or 'quadratic' recognised ephemeris types")

        self.coeff  = [float(s) for s in subv[2::2]]
        self.ecoeff = [float(s) for s in subv[3::2]]
        if self.time == 'HJD' or self.time == 'BJD':
            if self.coeff[0] < 2400000.:
                raise Exception('Ephem: for HJD or BJD expect T0 > 2400000')
        elif self.time == 'HMJD' or self.time == 'BMJD':
            if self.coeff[0] > 70000.:
                raise Exception('Ephem: for HMJD or BMJD expect T0 < 70000')
        else:
            raise Exception('Ephem: recognised times are HJD, BJD, HMJD or BMJD')

    def phase(self, time):
        """
        Computes phase corresponding to a given time
        """
        pnew = (time - self.coeff[0]) / self.coeff[1]
        if self.poly == 'quadratic':
            pold = pnew - 1
            while np.abs(pnew-pold) > 1.e-8:
                pold = pnew
                pnew = (time - self.coeff[0] - self.coeff[1]*pold**2) / \
                       self.coeff[1]
        return pnew

    def etime(self, cycle):
        """
        Computes uncertainty in time of ephemeris at a given phase
        """
        esum = 0
        fac  = 1.
        for ecf in self.ecoeff:
            esum += fac*ecf**2
            fac  *= cycle**2
        return m.sqrt(esum)


class Sdata (object):
    """
    Stores positional, ephemeris and line weight data on a star.
    """
    def __init__(self, ra, dec, ephem, lweight):
        self.ra  = ra
        self.dec = dec
        self.eph = ephem
        self.lw  = lweight

class Switch (object):
    """
    Stores switch target data.
    """
    def __init__(self, line):
        name, ut, delta = line.split('|')
        self.name = name.strip()
        utv = [int(i) for i in ut.split(':')]
        utc = float(utv[0])
        if len(utv) > 1:
            utc += float(utv[1])/60.
        if len(utv) > 2:
            utc += float(utv[2])/3600.
        self.utc   = utc
        self.delta = float(delta)/60.

class Prange(object):
    """
    Stores phase or time range data
    Very little to this class.
    """
    def __init__(self, name):
        self.name   = name
        self.prange = []

    def add(self, line):
        """
        Given a line of information with 4 components,
        a phase or time range, a pgpplot colour index and a line
        width, this stores the parameters in a an internal list
        prange. Distinguish phase from time (JD) by > or < 1000.
        """
        p1, p2, col, lw = line.split()
        p1  = float(p1)
        p2  = float(p2)
        if p1 < 1000.:
            p2  = p2 - m.floor(p2-p1)
            p_or_t = 'Phase'
        else:
            p_or_t = 'Time'

        col = int(col)
        lw  = int(lw)
        self.prange.append([p1, p2, col, lw, p_or_t])

def tnt_alert(alt, az):
    """
    TNT horizon is complicated by the TV mast. This warns
    that one is being obscured.

    Arguments::

      alt : altitude, degrees
      az  : azimuth, degrees

    Returns with True if the mast is in the way.
    """
    # three critical points defining mast are in (alt,az) (21.0,25.5) (lower
    # left-hand point), (73.5,33.5) (apex), (21.0,50.) (lower right-hand
    # point). Assume that these define two great circles, and find the axial
    # vectors of these great circles.
    ALT  = np.radians(np.array([21.0,73.5,21.0]))
    AZ   = np.radians(np.array([25.5,33.5,50.0]))
    calt, salt = np.cos(ALT), np.sin(ALT)
    caz,  saz  = np.cos(AZ), np.sin(AZ)
    v1  = subs.Vec3(saz[0]*calt[0], caz[0]*calt[0], salt[0])
    v2  = subs.Vec3(saz[1]*calt[1], caz[1]*calt[1], salt[1])
    v3  = subs.Vec3(saz[2]*calt[2], caz[2]*calt[2], salt[2])

    # a1, a2 defined to be downward pointing axial vectors corresponding to
    # great cirles representing each extreme of the mast. Inside mast if
    # actual vector gives positive dot product with both axial vectors.
    a1  = subs.cross(v1,v2)
    a2  = subs.cross(v2,v3)

    ralt, raz = m.radians(alt), m.radians(az)
    v = subs.Vec3(m.sin(raz)*m.cos(ralt), m.cos(raz)*m.cos(ralt), m.sin(ralt))

    return subs.dot(a1,v) > 0 and subs.dot(a2,v) > 0

