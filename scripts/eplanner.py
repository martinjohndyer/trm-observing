#!/usr/bin/env python

from __future__ import print_function

usage = \
"""
Plots out observing schedules for eclipsers. You supply a file of target
positions and ephemerides, another defining the phase ranges of interest and
optionally a third specifying when to switch targets, and this will make a
graphical representation of the results. It also allows you to plot a line
representing the objects you wish to observe.

The tracks for the stars start and stop when the Sun is 10 degrees below the
horizon, although that is really pushing things. Don't expect to be able to
push beyond these limits or even close to them in most cases.

If you see horizontal black error bars at the far right, these indicate +/- 1
sigma uncertainties on ephemerides. As quoted uncertainties can barely ever be
trusted, they are indicative only. If they are big, be careful: your eclipse
may not appear when you expect.

Grey boxes represent zenith holes for Alt/Az telescopes. Some like the VLT
can't track close to the zenith, and you may not be able to observe during
these intervals.
"""
import argparse
import datetime
import math as m

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy import time, coordinates as coord, units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz

from trm import observing

# Longitude, latitude, height keyed by name
SITES = {
    'WHT' : ('-17 52 53.9', '+28 45 38.3', 2332.),
    'NTT' : ('-70 44 00.', '-29 15 00.', 2400.),
    'VLT' : ('-70 24 9.9', '-24 37 30.3', 2635.),
    'TNT' : ('+98 28 00', '+18 34 00',2457.),
    }


if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(
        description=usage, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional
    parser.add_argument(
        'stardata',
        help='file containing positions and ephemerides. It should have the same format of my C++ ephemeris file.')

    parser.add_argument(
        'telescope', help='Telescope name, e.g. WHT, NTT, VLT, TNT')

    parser.add_argument(
        'pranges',
        help='file of stars to include and phase ranges of interest. Each entry should start with the name of the star on its own on a line, and then follow it with the phase ranges in the form:\np1 p2 col lw\nwhere p1 p2 is the range of phase to indicate, col is the pgplot colour index (2 = red, etc), lw is the pgplot line width. It is OK to have p1 > p2; the programme will map it into a positive range.')

    parser.add_argument(
        'date', help='date at start of night, format: "YYYY-MM-DD"')

    # optional
    parser.add_argument(
        '-d', dest='device', default='/xs',
        help='PGPLOT plot device, e.g. plot.ps/cps, plot.ps/vcps.')

    parser.add_argument(
        '-t', dest='twilight', type=float, default=-15,
        help='how many degrees below horizon for twilight.')

    parser.add_argument(
        '-c', dest='csize', type=float, default=1.3,
        help='default character size for star names.')

    parser.add_argument(
        '-a', dest='airmass', type=float, default=2,
        help='airmass limit.')

    parser.add_argument(
        '-f', dest='fudge', type=float, default=0.55,
        help='fudge factor to allow adjustment of space at left for star names.')

    parser.add_argument(
        '-o', dest='offset', type=float, default=0.2,
        help='offset as fraction of plot width at which to print duplicate names.')

    parser.add_argument(
        '-p', dest='prefix', type=argparse.FileType('r'),
        default=None,
        help='target name prefixes. These are placed before the left-hand list of target names to allow one to add e.g. priorities.')

    parser.add_argument(
        '-s', dest='switch', type=argparse.FileType('r'), default=None,
        help='switch targets data file. Each line should have the format:\nstar name | start switch time (e.g. 12:34) | dead time\nwhere the dead time is the time taken (in minutes) to make the switch. The line will be split on the pipe | characters which must therefore not appear in the star name. Use a star name = None as the final switch to terminate early. The times should increase monotonically, so use times > 24 if necessary.')

    # parse them
    args = parser.parse_args()

    # Interpret date
    date = time.Time(args.date, out_subfmt='date')

    # Location
    if args.telescope in SITES:
        info = SITES[args.telescope]
        site = coord.EarthLocation.from_geodetic(info[0],info[1],info[2])
    else:
        site = coord.EarthLocation(args.telescope)

    # Load position and ephemeris data
    peinfo = observing.load_pos_eph(args.stardata)

    # Load phase / time ranges
    prinfo = observing.load_ptranges(args.pranges, peinfo)

    # Load switches
    swinfo = observing.load_switches(args.switch, peinfo)

    # Load prefixes
    prefixes = observing.load_prefixes(args.prefix, peinfo)

    # Compute maximum length of the names
    if len(prefixes):
        mpre = max([len(n) for n in prefixes.values()]) + 1
    else:
        mpre = 0

    left = max([len(n) for n in peinfo.keys()]) + mpre

    # Rather primitive times to define sun down and up; won't work
    # properly in far north or south or around the dateline, but should be
    # OK for most sites.

    # 'date' is set at the UTC of the start of the entered date, i.e. 0 hours. To get to the
    # day time before the night of interest, we advance by 0.5 days and apply an
    # offset in longitude. This should be the local mid-day.
    toffset = time.TimeDelta(0.5 - site.longitude.degree/360., format='jd')

    day1  = date + toffset
    night = day1 + time.TimeDelta(0.5, format='jd')
    day2  = night + time.TimeDelta(0.5, format='jd')

    # sunset & sunrise times
    sunset = observing.sun_at_alt(day1, night, site, -0.25)
    sunrise = observing.sun_at_alt(night, day2, site, -0.25)
    sunset.out_subfmt='date_hm'
    sunrise.out_subfmt='date_hm'
    print('Sun is at alt=-0.25 at {0:s} and {1:s}'.format(sunset.iso,sunrise.iso))

    # Sun should at least be at -10 really
    rset = observing.sun_at_alt(day1, night, site, -10.)
    rrise = observing.sun_at_alt(night, day2, site, -10.)
    rset.out_subfmt='date_hm'
    rrise.out_subfmt='date_hm'
    print('Sun is at alt=-10.0 at {0:s} and {1:s}'.format(rset.iso,rrise.iso))

    # User-defined twilight
    twiset = observing.sun_at_alt(day1, night, site, args.twilight)
    twirise = observing.sun_at_alt(night, day2, site, args.twilight)
    twiset.out_subfmt='date_hm'
    twirise.out_subfmt='date_hm'
    print('Sun is at alt={0:5.1f} at {1:s} and {2:s}'.format(args.twilight,twiset.iso,twirise.iso))

    # integer offset
    isun = int(sunset.mjd)

    cols = {
        2 : (0.7,0,0),
        3 : (0,0.7,0),
        4 : (0,0,0.7),
        5 : (0.85,0.85,0.85),
        9 : (0.8,0.8,0.8),
    }

    # Start the plot
    fig = plt.figure()
    ax = plt.axes([0.2,0.1,0.75,0.85])

    # sunrise / sunset
    utc1 = 24.*(sunset.mjd-isun)
    utc2 = 24.*(sunrise.mjd-isun)

    # correct times to lie within range
    for sw in swinfo:
        if sw.utc < utc1 and sw.utc+24 < utc2:
            sw.utc += 24.

    # rise / set as defined by user-defined altitude
    utc3 = 24.*(twiset.mjd-isun)
    utc4 = 24.*(twirise.mjd-isun)

    # rise / set when at -10 which is the absolute limit I think.
    utc5 = 24.*(rset.mjd-isun)
    utc6 = 24.*(rrise.mjd-isun)

    # set up general scale, draw vertical dashed lines every hour
    utstart = utc1 - 0.4
    utend = utc2 + 0.4

    n1 = int(np.ceil(utc1))
    n2 = int(np.ceil(utc2))

    for n in range(n1,n2):
        plt.plot([n,n],[0,1],'--',color=cols[5])

    # mark sunset / twiset / twirise / sunrise
    kwargs = {'color' : cols[2], 'horizontalalignment' : 'center'}
    plt.text(utc1, 1.02, 'sunset', **kwargs)
    plt.text(utc2, 1.02, 'sunrise', **kwargs)
    plt.text(utc3, 1.02, str(args.twilight), **kwargs)
    plt.text(utc4, 1.02, str(args.twilight), **kwargs)

    kwargs = {'color' : cols[2]}
    plt.plot([utc1,utc1],[0,1],'--',**kwargs)
    plt.plot([utc2,utc2],[0,1],'--',**kwargs)
    plt.plot([utc3,utc3],[0,1],'--',**kwargs)
    plt.plot([utc4,utc4],[0,1],'--',**kwargs)

    # 400 points from start to end
    utcs = np.linspace(utc5,utc6,400)
    mjds = isun + utcs/24.
    mjd6 = isun+utc6/24.

    # Loop through the stars listed in the ephemeris file, computing their hour
    # angles at midnight.  This establishes the plot order.
    has = []
    midnight = time.Time((sunset.mjd+sunrise.mjd)/2.,format='mjd')
#   for key, star in peinfo.items():
#       airmasses, alt, az, ha, pa, delz = \
#           sla.amass(midnight,longit,latit,height,star.ra,star.dec)
#       has.append((ha, key))
#
#    keys = [item[1] for item in sorted(has,key=lambda ha: ha[0], reverse=True)]

    # Now compute plot positions
#    ys = {}
#    for i,key in enumerate(keys):
#        ys[key] = (len(keys)-i)/float(len(keys)+1)


    plt.xlabel('UTC')
    print(date,args.telescope)
    plt.title('Night starting {0!s} at the {1:s}, airmass < {2:3.1f}'.format(
            date, args.telescope, args.airmass))

    plt.xlim(utstart, utend)
    plt.ylim(0,1.05)

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().set_visible(False)

    # Modify hours labels to lie in range 1 to 24
    fig.canvas.draw()
    labels = [item.get_text() for item in ax.get_xticklabels()]
    for i in range(len(labels)):
        try:
            if int(labels[i]) > 24:
                labels[i] = str(int(labels[i])-24)
        except:
            pass
    ax.set_xticklabels(labels)

    plt.show()
    exit(1)



#    pgsvp(xv1, xv2, yv1+ch/40, yv2-ch/40)
#    pgswin(utstart,utend,0,1)



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
        airmasses, alts, azs, has, pa, delz = \
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

            ok = (mjds > mjd_start) & (mjds < mjd_end)

            if tel == 'TNT':
                # TNT specific, first for the aerial
                start = True
                aerial = []
                for alt,az,utc in zip(alts[ok],azs[ok],utcs[ok]):
                    if az < 0.: az += 360.

                    if observing.tnt_alert(alt, az):
                        if start:
                            # start bad period
                            air_start = utc
                            air_end   = utc
                            start     = False
                        else:
                            # update end time of bad period
                            air_end = utc
                    else:
                        if not start:
                            # We are out of it. Record bad period.
                            aerial.append((air_start,air_end))
                            start = True

                if not start:
                    # We were still in a bad period at the end. Record it.
                    aerial.append((air_start,air_end))

                # Plots the bad periods
                for air_start, air_end in aerial:
                    pgsci(9)
                    pgsfs(1)
                    pgrect(air_start,air_end,y-0.01,y+0.01)
                    pgsfs(2)
                    pgsls(1)
                    pgsci(1)
                    pgrect(air_start,air_end,y-0.01,y+0.01)

                # then for close to zenith
                start = True
                end   = False
                air_start, air_end = utc_end, utc_start
                for alt, utc in zip(alts[ok], utcs[ok]):
                    if start and alt > 88.:
                        air_start = utc
                        air_end   = utc_end
                        start     = False

                    if not start and not end and alt < 88.:
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

            elif tel == 'WHT' or tel == 'NTT':
                # WHT / NTT specific
                start = True
                end   = False
                air_start, air_end = utc_end, utc_start
                for alt, utc in zip(alts[ok], utcs[ok]):
                    if start and alt > 87.:
                        air_start = utc
                        air_end   = utc_end
                        start     = False

                    if not start and not end and alt < 87.:
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

            if key in prinfo:
                # handles the time ranges only
                pr   = prinfo[key]

                pgsls(1)
                pranges = pr.prange
                for p1, p2, col, lw, p_or_t in pranges:
                    if p_or_t == 'Time':
                        utc1, utc2 = 24.*(p1-isun), 24.*(p2-isun)
                        if utc1 < utc_end and utc2 > utc_start:
                            ut1  = max(utc1, utc_start)
                            ut2  = min(utc2, utc_end)
                            if ut1 < ut2:
                                pgsci(col)
                                pgslw(lw)
                                pgmove(ut1, y)
                                pgdraw(ut2, y)

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
                    for p1, p2, col, lw, p_or_t in pranges:
                        if p_or_t == 'Phase':
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

                # draws dots at phase zero
                d1 = np.ceil(pstart)
                d2 = np.floor(pend)
                nphs = int(np.ceil(d2 - d1))+1
                for n in range(nphs):
                    ut = utc_start + (utc_end-utc_start)*(d1 + n - pstart)/(pend-pstart)
                    pgsci(1)
                    pgslw(3)
                    pgsch(1)
                    pgpt1(ut, y, 17)

            # draw vertical bar at meridian
            if ok.any():
                hamin, hamax = has[ok].min(), has[ok].max()
                if hamin < 0 and hamax > 0.:
                    has = np.abs(has[ok])
                    utc_mer = utcs[ok][has.argmin()]
                    pgslw(4)
                    pgsci(1)
                    pgmove(utc_mer, y-1.3*lbar)
                    pgdraw(utc_mer, y+1.3*lbar)
                    pgslw(1)

    # add target names at left
    pgslw(1)
    pgsci(1)
    pgsvp(0,1,yv1+ch/40,yv2-ch/40)
    pgswin(0,1,0,1)
    pgsch(args.csize)
    for i,key in enumerate(keys):
        y = (len(peinfo)-i)/float(len(peinfo)+1)
        if key in prefixes:
            name = prefixes[key] + ' ' + key
        elif len(prefixes):
            name = mpre*' ' + key
        else:
            name = key

        pgslw(peinfo[key].lw)
        pgptxt(0.,(len(peinfo)-i)/float(len(peinfo)+1),0,0,name)

    pgclos()

    print('Airmass limit  =',args.airmass)
    print('Twilight limit =',args.twilight,'degrees')
