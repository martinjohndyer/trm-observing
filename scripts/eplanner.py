#!/usr/bin/env python

usage = \
"""
Plots out observing schedule for eclipsers. You supply a file of target
positions and ephemerides, another defining the phase ranges of interest and
optionally a third specifying wqhen to switch targets, and this will make a
graphical representation of the results. It also allows you to plot a line
representing the objects you wish to observe.

The tracks for the stars start and stop when the Sun is 10 degrees below the
horizon, although that is really pushing things. Don't expect to be able to
push beyond these limits.

If you see horizontal black error bars at the far right, these indicate +/- 1
sigma uncertainties on ephemerides. As quoted uncertainties can barely ever be
trusted, they are indicative only. If they are big, be careful: your eclipse
may not appear when you expect.
"""

import math as m
import numpy as np
import argparse
from ppgplot import *
from trm import subs, sla, observing

# arguments
parser = argparse.ArgumentParser(description=usage,\
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# positional
parser.add_argument('stardata', type=argparse.FileType('r'),
                   help='general file containing positions and ephemerides. It should have the same format of my C++ ephemeris file.')

parser.add_argument('telescope', help='Telescope name, e.g. WHT, VLT')

parser.add_argument('pranges', type=argparse.FileType('r'),
                    help='file of stars to include and phase ranges of interest. Each entry should start with the name of the star on its own on a line, and then follow it with the phase ranges in the form:\np1 p2 col lw\nwhere p1 p2 is the range of phase to indicate, col is the pgplot colour index (2 = red, etc), lw is the pgplot line width. It is OK to have p1 > p2; the programme will map it into a positive range.')

parser.add_argument('date', help='date at start of night, e.g. "01 Aug 2012"')

# optional
parser.add_argument('-d', dest='device', default='/xs',
                   help='PGPLOT plot device, e.g. plot.ps/cps, plot.ps/vcps.')

parser.add_argument('-t', dest='twilight', type=float, default=-15,
                   help='how many degrees below horizon for twilight.')

parser.add_argument('-c', dest='csize', type=float, default=1.3,
                   help='default character size for star names.')

parser.add_argument('-a', dest='airmass', type=float, default=2,
                   help='airmass limit.')

parser.add_argument('-f', dest='fudge', type=float, default=0.55,
                   help='fudge factor to allow adjustment of space at left for star names.')

parser.add_argument('-o', dest='offset', type=float, default=0.2,
                   help='offset as fraction of plot width at which to print duplicate names.')

parser.add_argument('-s', dest='switch', type=argparse.FileType('r'), default=None,
                   help='switch targets data file. Each line should have the format:\nstar name | start switch time (e.g. 12:34) | dead time\nwhere the dead time is the time taken (in minutes) to make the switch. The line will be split on the pipe | characters which must therefore not appear in the star name. Use a star name = None as the final switch to terminate early. The times should increase monotonically, so use times > 24 if necessary.')

# parse them
args = parser.parse_args()

# Interpret date
year,month,day = subs.date2dmy(args.date)

# Get observatory parameters
tel,obs,longit,latit,height = subs.observatory(args.telescope)

# Load position and ephemeris data
peinfo = {}
count = 0
nline = 0
name  = None
for line in args.stardata:
    nline += 1
    try:
        if not line.startswith('#') and not line.isspace():
            count += 1
            if count == 1:
                name = line.strip()
            elif count == 2:
                ra,dec,system = subs.str2radec(line.strip())
            elif count == 3:
                try:
                    eph = observing.Ephemeris(line)
                    peinfo[name] = observing.Sdata(ra, dec, eph)
                except:
                    print 'No valid ephemeris data found for',name
                    peinfo[name] = Sdata(ra, dec, None)
                count = 0
    except Exception, err:
        print err
        print 'Line number',nline
        print 'Line =',line.strip()
        if name:
            print 'Name = ' + name
        else:
            print 'Name undefined'
        print 'Program aborted.'
        exit(1)

args.stardata.close()
print 'Data on',len(peinfo),'stars loaded.'


# Load phase ranges
prinfo = {}
count = 0
for line in args.pranges:
    try:
        if line.startswith('#') or line.isspace():
            if count:
                if name in peinfo:
                    prinfo[name] = pr
                else:
                    print name,'not found in position/ephemeris file and will be skipped.'
            count = 0
        else:
            count += 1
            if count == 1:
                name = line.strip()
                pr   = observing.Prange(name)
            elif count > 1:
                p1, p2, col, lw = line.split()
                p1  = float(p1)
                p2  = float(p2)
                p2  = p2 - m.floor(p2-p1)
                col = int(col)
                lw  = int(lw)
                pr.prange.append([p1, p2, col, lw])
    except ValueError:
        print 'ValueError found in',line
        exit(1)
args.pranges.close()
print 'Data on',len(prinfo),'phase ranges loaded.'

# Load switches
swinfo = []
if args.switch is not None:
    first = True
    for line in args.switch:
        if not line.startswith('#') and not line.isspace():
            swinfo.append(observing.Switch(line))
            if swinfo[-1].name != 'None' and swinfo[-1].name not in peinfo:
                raise Exception('switch star: ' + swinfo[-1].name + \
                                ' not found in position/ephemeris file.')
            if first:
                utold = swinfo[-1].utc
            else:
                utnew = swinfo[-1].utc
                if utnew < utold:
                    raise Exception('switch: times not increasing.')
                utold = utnew
    args.switch.close()
    print len(swinfo),'target switches loaded.'
else:
    print 'No target switches loaded.'

# Estimate the maximum length of the names
left = max([len(n) for n in prinfo.keys()])

# Rather primitive times to define sun down and up; won't work
# properly in far north or south or around the dateline
mjd_day1  = sla.cldj(year, month, day) - longit/360. + 0.5
mjd_night = mjd_day1 + 0.5
mjd_day2  = mjd_day1 + 1.0

sunset  = sla.sun_at_elev(longit, latit, height, mjd_day1, mjd_night, -0.25)
sunrise = sla.sun_at_elev(longit, latit, height, mjd_night, mjd_day2, -0.25)

twiset  = sla.sun_at_elev(longit, latit, height, mjd_day1, mjd_night,
                          args.twilight)
twirise = sla.sun_at_elev(longit, latit, height, mjd_night, mjd_day2,
                          args.twilight)

# really stop at these points
rset  = sla.sun_at_elev(longit, latit, height, mjd_day1, mjd_night, -10.)
rrise = sla.sun_at_elev(longit, latit, height, mjd_night, mjd_day2, -10.)

# integer offset
isun = int(sunset)

# Start the plot
pgopen(args.device)

# set up generalities
pgslw(3)
pgscr(0,1,1,1)
pgscr(1,0,0,0)
pgscr(2,0.7,0,0)
pgscr(3,0,0.7,0)
pgscr(4,0,0,0.7)
pgscr(9,0.8,0.8,0.8)
pgsch(1.3)
ch = pgqch()

pgsvp(0,1,0,1)
xr1,xr2,yr1,yr2 = pgqvp(2)
aspect = (yr2-yr1)/(xr2-xr1)

# compute viewport
xv1  = args.fudge*aspect*left*args.csize/40.
xv2  = 1.-aspect*2.*ch/40.
yv1  = 4.*ch/40.
yv2  = 1.-4.*ch/40.

# sunrise / sunset
utc1 = 24.*(sunset-isun)
utc2 = 24.*(sunrise-isun)

# correct times to lie within range
for sw in swinfo:
    if sw.utc < utc1 and sw.utc+24 < utc2:
        sw.utc += 24.

# rise / set as defined by user-defined angle
utc3 = 24.*(twiset-isun)
utc4 = 24.*(twirise-isun)

# rise / set when at -12 which is really the limit I think.
utc5 = 24.*(rset-isun)
utc6 = 24.*(rrise-isun)

# set up general scale, draw vertical dashed lines every hour
utstart = utc1
utend   = utc2 + 0.4
pgsvp(xv1,xv2,yv1,yv2)
pgswin(utstart,utend,0,1)
n1 = int(np.ceil(utc1))
n2 = int(np.ceil(utc2))
pgsls(2)
pgscr(5,0.85,0.85,0.85)
pgsci(5)
for n in range(n1,n2):
    pgmove(n,0)
    pgdraw(n,1)

# draw axes (bit fiddly to allow 24 to 0 transition)
pgsls(1)
if utc2 > 24:
    pgsci(1)
    pgsvp(xv1,xv1+(xv2-xv1)*(24.-utstart)/(utend-utstart),yv1,yv2)
    pgswin(utstart,23.9999,0,1)
    pgsci(4)
    pgbox('bcnst',1,2,'',0,0)
    pgsvp(xv1+(xv2-xv1)*(24.-utstart)/(utend-utstart),xv2,yv1,yv2)
    pgswin(0,utend-24,0,1)
    pgbox('bcnst',1,2,'',0,0)

    # back to general scale
    pgsvp(xv1,xv2,yv1,yv2)
    pgswin(utstart,utend,0,1)
else:
    pgsci(1)
    pgbox('bcnst',1,0,'',0,0)

# mark sunset / twiset / twirise / sunrise
pgsci(2)
pgsch(0.9)
pgptxt(utc1, 1.02, 0, 0.5, 'sunset')
pgptxt(utc2, 1.02, 0, 0.5, 'sunrise')
pgptxt(utc3, 1.02, 0, 0.5, str(args.twilight))
pgptxt(utc4, 1.02, 0, 0.5, str(args.twilight))
pgsls(2)
pgmove(utc1,0)
pgdraw(utc1,1)
pgmove(utc2,0)
pgdraw(utc2,1)
pgmove(utc3,0)
pgdraw(utc3,1)
pgmove(utc4,0)
pgdraw(utc4,1)
pgsch(1.5)
pglab('UTC',' ',args.date + ' (' + args.telescope + ', airmass < ' + \
      str(args.airmass) + ')')

utcs = np.linspace(utc5,utc6,400)
mjds = isun + utcs/24.
mjd6 = isun+utc6/24.

pgsvp(xv1, xv2, yv1+ch/40, yv2-ch/40)
pgswin(utstart,utend,0,1)

# Loop through the stars listed in the ephemeris file, computing their hour
# angles at midnight.  This establishes the plot order.
has = []
midnight = (sunset+sunrise)/2.
for key, star in peinfo.iteritems():
    airmasses, alt, az, ha, pa, delz = \
        sla.amass(midnight,longit,latit,height,star.ra,star.dec)
    has.append((ha, key))

keys = [item[1] for item in sorted(has,key=lambda ha: ha[0], reverse=True)]

# Now compute plot positions
ys = {}
for i,key in enumerate(keys):
    ys[key] = (len(keys)-i)/float(len(keys)+1)

# Plot switches
if args.switch is not None:
    pgsls(1)
    pgscr(5,0.8,0.8,0.8)
    pgslw(10)
    pgsci(5)

    first = True
    for sw in swinfo:
        if first:
            yold = ys[sw.name]
            pgmove(sw.utc, yold)
            first = False
        else:
            if sw.name == 'None':
                pgdraw(sw.utc, yold)
                break
            else:
                pgdraw(sw.utc, yold)
                yold = ys[sw.name]
                pgdraw(sw.utc+sw.delta, yold)

# Loop through the stars listed in the phase ranges.
for i, key in enumerate(keys):
    star = peinfo[key]
    y    = ys[key]

    # Compute airmasses etc for all points
    airmasses, alts, azs, ha, pa, delz = \
        sla.amass(mjds,longit,latit,height,star.ra,star.dec)

    start = True
    end   = False
    utc_start, utc_end = utc6, utc5
    for airmass,utc,mjd in zip(airmasses,utcs,mjds):

        if start and airmass < args.airmass:
            utc_start = utc
            utc_end   = utc6
            mjd_start = mjd
            mjd_end   = mjd6
            start     = False

        if not start and not end and airmass > args.airmass:
            utc_end = utc
            mjd_end = mjd
            end     = True
            break

    pgsls(2)
    pgsci(1)
    pgslw(1)
    if utc_start < utc_end:
        # Draw line from rise to set airmass
        pgmove(utc_start, y)
        pgdraw(utc_end, y)
        pgsls(1)
        lbar = min(0.01, 1./len(keys)/4.)
        pgmove(utc_start, y-lbar)
        pgdraw(utc_start, y+lbar)
        pgmove(utc_end, y-lbar)
        pgdraw(utc_end, y+lbar)

        if utc_start > utstart + args.offset*(utend-utstart):
            # repeat target name if the start is late
            pgscf(1)
            pgsch(0.7*args.csize)
            pgptxt(utc_start-0.2,y,0,1.0,key)

        if tel == 'TNT':
            # TNT specific
            ok = (mjds > mjd_start) & (mjds < mjd_end)
            start = True
            end   = False
            air_start, air_end = utc_end, utc_start
            for alt,az,utc in zip(alts[ok],azs[ok],utcs[ok]):
                if az < 0.: az += 360.

                if start and observing.tnt_alert(alt, az):
                    air_start = utc
                    air_end   = utc_end
                    start     = False

                if not start and not end and not observing.tnt_alert(alt, az):
                    air_end = utc
                    end     = True
                    break

            if air_start < air_end:
                pgsci(9)
                pgsfs(1)
                pgrect(air_start,air_end,y-0.01,y+0.01)
                pgsfs(2)
                pgsls(1)
                pgsci(1)
                pgrect(air_start,air_end,y-0.01,y+0.01)

        elif tel == 'VLT':
            # VLT specific
            ok = (mjds > mjd_start) & (mjds < mjd_end)
            start = True
            end   = False
            air_start, air_end = utc_end, utc_start
            for alt, utc in zip(alts[ok], utcs[ok]):
                if start and alt > 86.:
                    air_start = utc
                    air_end   = utc_end
                    start     = False

                if not start and not end and alt < 86.:
                    air_end = utc
                    end     = True
                    break

            if air_start < air_end:
                pgsci(9)
                pgsfs(1)
                pgrect(air_start,air_end,y-0.01,y+0.01)
                pgsfs(2)
                pgsls(1)
                pgsci(1)
                pgrect(air_start,air_end,y-0.01,y+0.01)


        # Compute phase info
        tt,tdb,btdb_start,hutc_start,htdb,vhel,vbar = \
            sla.utc2tdb(mjd_start,longit,latit,height,star.ra,star.dec)

        tt,tdb,btdb_end,hutc_end,htdb,vhel,vbar = \
            sla.utc2tdb(mjd_end,longit,latit,height,star.ra,star.dec)

        if star.eph:
            eph = star.eph

            if eph.time == 'HJD':
                pstart = eph.phase(hutc_start + 2400000.5)
                pend   = eph.phase(hutc_end + 2400000.5)
            elif eph.time == 'HMJD':
                pstart = eph.phase(hutc_start)
                pend   = eph.phase(hutc_end)
            elif eph.time == 'BMJD':
                pstart = eph.phase(btdb_start)
                pend   = eph.phase(btdb_end)
            elif eph.time == 'BJD':
                pstart = eph.phase(btdb_start + 2400000.5)
                pend   = eph.phase(btdb_end + 2400000.5)
            else:
                raise Exception('Unrecognised type of time = ' + eph.time)

            # Compute uncertainty in predictions
            delta = 24.*min(eph.etime((pstart+pend)/2.), eph.coeff[1]/2.)
            if delta > 0.03:
                pgsls(1)
                pgslw(5)
                pgsci(1)
                pgmove(utend-2.*delta,y)
                pgdraw(utend,y)
                pgmove(utend-2.*delta,y-0.005)
                pgdraw(utend-2.*delta,y+0.005)
                pgmove(utend,y-0.005)
                pgdraw(utend,y+0.005)

            if key in prinfo:
                pr   = prinfo[key]

                # Draw phase ranges of interest
                pgsls(1)
                pranges = pr.prange
                for p1, p2, col, lw in pranges:
                    d1    = pstart + (p1 - pstart) % 1 - 1
                    d2    = pend   + (p1 - pend) % 1
                    nphs  = int(np.ceil(d2 - d1))
                    for n in range(nphs):
                        ut1  = utc_start + (utc_end-utc_start)*\
                               (d1 + n - pstart)/(pend-pstart)
                        ut2  = ut1 + (utc_end-utc_start)/(pend-pstart)*(p2-p1)
                        ut1  = max(ut1, utc_start)
                        ut2  = min(ut2, utc_end)
                        if ut1 < ut2:
                            pgsci(col)
                            pgslw(lw)
                            pgmove(ut1, y)
                            pgdraw(ut2, y)

# add target names at left
pgslw(1)
pgsci(1)
pgsvp(0,1,yv1+ch/40,yv2-ch/40)
pgswin(0,1,0,1)
pgsch(args.csize)
for i,key in enumerate(keys):
    y = (len(peinfo)-i)/float(len(peinfo)+1)
    pgptxt(0.,(len(peinfo)-i)/float(len(peinfo)+1),0,0,key)

pgclos()

print 'Airmass limit  =',args.airmass
print 'Twilight limit =',args.twilight,'degrees'
