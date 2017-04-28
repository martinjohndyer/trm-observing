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
        '-f', dest='hcopy', default=None,
        help='Name for hard copy file, e.g. plot.pdf, plot.png (there will be no interactive plot in this case)')

    parser.add_argument(
        '-t', dest='twilight', type=float, default=-15,
        help='how many degrees below horizon for twilight.')

    parser.add_argument(
        '-c', dest='csize', type=float, default=1.0,
        help='default character size for star names.')

    parser.add_argument(
        '-a', dest='airmass', type=float, default=2,
        help='airmass limit.')

    parser.add_argument(
        '-o', dest='offset', type=float, default=0.2,
        help='offset as fraction of plot width at which to print duplicate names.')

    parser.add_argument(
        '-d', dest='divide', type=float, default=0.2,
        help='divide point between names on left and plot on right as fraction of total width')

    parser.add_argument(
        '--hfrac', type=float, default=0.8,
        help='plot height as fraction of total')

    parser.add_argument(
        '--bfrac', type=float, default=0.1,
        help='offset of bottom of plot as fraction of total')

    parser.add_argument(
        '-p', dest='prefix', default=None,
        help='target name prefixes. These are placed before the left-hand list of target names to allow one to add e.g. priorities.')

    parser.add_argument(
        '-s', dest='switch', default=None,
        help='switch targets data file name. Each line should have the format:\nstar name | start switch time (e.g. 12:34) | dead time\nwhere the dead time is the time taken (in minutes) to make the switch. The line will be split on the pipe | characters which must therefore not appear in the star name. Use a star name = None as the final switch to terminate early. The times should increase monotonically, so use times > 24 if necessary.')

    parser.add_argument(
        '-x', dest='width', type=float, default=8,
        help='plot width, inches')

    parser.add_argument(
        '-y', dest='height', type=float, default=5,
        help='plot height, inches')

    # parse them
    args = parser.parse_args()

    assert(args.airmass > 1 and args.airmass < 6)
    assert(args.width > 0 and args.height > 0)

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

    # Rather primitive times to define Sun down and up; won't work
    # properly in far north or south or around the dateline, but should be
    # OK for most sites.

    # 'date' is set at the UTC of the start of the entered date, i.e. 0
    # hours. To get to the day time before the night of interest, we advance
    # by 0.5 days and apply an offset in longitude. This should be the local
    # mid-day, then half day steps are used to locate midnight and the
    # following mid day. Armed with these we can narrow down on the sunset /
    # rise etc times.
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

    # the lines in the plot will terminate at the point the Sun is at -10 as the absolute
    # limit for observation.
    rset = observing.sun_at_alt(day1, night, site, -10.)
    rrise = observing.sun_at_alt(night, day2, site, -10.)
    rset.out_subfmt='date_hm'
    rrise.out_subfmt='date_hm'
    print('Sun is at alt=-10.0 at {0:s} and {1:s}'.format(rset.iso,rrise.iso))

    # user-defined twilight -- will be indicated on the plot.
    twiset = observing.sun_at_alt(day1, night, site, args.twilight)
    twirise = observing.sun_at_alt(night, day2, site, args.twilight)
    twiset.out_subfmt='date_hm'
    twirise.out_subfmt='date_hm'
    print('Sun is at alt={0:5.1f} at {1:s} and {2:s}'.format(args.twilight,twiset.iso,twirise.iso))

    # integer offset
    isun = int(sunset.mjd)

    # simulating the old PGPLOT colours
    cols = {
        2 : (0.7,0,0),
        3 : (0,0.7,0),
        4 : (0,0,0.7),
        5 : (0.85,0.85,0.85),
        9 : (0.8,0.8,0.8),
    }

    # Start the plot
    fig = plt.figure(figsize=(args.width,args.height))

    # left and right axes. Left for labels, right for the main stuff
    axl = plt.axes([0.01, args.bfrac, args.divide, args.hfrac])
    axr = plt.axes([args.divide, args.bfrac, 0.99-args.divide, args.hfrac])

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
        axr.plot([n,n],[0,1],'--',color=cols[5])

    # mark sunset / twiset / twirise / sunrise
    kwargs = {'color' : cols[2], 'horizontalalignment' : 'center'}
    axr.text(utc1, 1.02, 'sunset', **kwargs)
    axr.text(utc2, 1.02, 'sunrise', **kwargs)
    axr.text(utc3, 1.02, str(args.twilight), **kwargs)
    axr.text(utc4, 1.02, str(args.twilight), **kwargs)

    kwargs = {'color' : cols[2]}
    axr.plot([utc1,utc1],[0,1],'--',**kwargs)
    axr.plot([utc2,utc2],[0,1],'--',**kwargs)
    axr.plot([utc3,utc3],[0,1],'--',**kwargs)
    axr.plot([utc4,utc4],[0,1],'--',**kwargs)

    # 400 points from start to end
    mjd6 = isun+utc6/24.
    utcs = np.linspace(utc5,utc6,400)
    mjds = isun + utcs/24.
    mjds = time.Time(mjds, format='mjd')

    # Loop through the stars listed in the ephemeris file, computing their hour
    # angles at midnight to establish the plot order
    has = []
    midnight = time.Time((sunset.mjd+sunrise.mjd)/2.,format='mjd')
    lst = midnight.sidereal_time(kind='mean', longitude=site.longitude)

    for key, star in peinfo.items():
        ha = (lst - star.position.ra).value
        if ha > 12.: ha -= 24
        if ha < -12.: ha += 24
        has.append((ha, key))

    keys = [item[1] for item in sorted(has, key=lambda ha: ha[0], reverse=True)]

    # Now compute plot positions
    ys = {}
    for i,key in enumerate(keys):
        ys[key] = (len(keys)-i)/float(len(keys)+1)


    # Plot switches
    if args.switch is not None:
        kwargs = {'color' : cols[5], 'lw' : 5}

        first = True
        for sw in swinfo:
            if first:
                xstart, ystart = sw.utc,ys[sw.name]
                first = False
            else:
                if sw.name == 'None':
                    plt.plot([xstart,sw.utc],[ystart,ystart],**kwargs)
                    break
                else:
                    xend, yend = sw.utc+sw.delta, ys[sw.name]
                    plt.plot([xstart, sw.utc, xend],[ystart, ystart, yend],**kwargs)
                    xstart, ystart = xend, yend


    # Loop through the stars listed in the phase ranges.
    for key in keys:
        star = peinfo[key]
        y = ys[key]

        # Compute airmasses for all times
        altazframe = AltAz(obstime=mjds, location=site)
        altaz = star.position.transform_to(altazframe)
        airmasses = altaz.secz.value
        alts = altaz.alt.value
        azs = altaz.az.value
        ok = alts > 90.-np.degrees(np.arccos(1./args.airmass))

        first = True
        afirst = True
        for flag, utc, mjd in zip(ok, utcs, mjds):
            if first and flag:
                ut_start = utc
                first = False
                if afirst:
                    mjd_first = mjd
                    utc_first = utc
                    afirst = False
            elif not flag and not first:
                first = True
                plt.plot([ut_start,utc],[y,y],'k--')
                lbar = min(0.01, 1./len(keys)/4.)
                plt.plot([ut_start,ut_start],[y-lbar,y+lbar],'k')
                plt.plot([utc,utc],[y-lbar,y+lbar],'k')
                mjd_last = mjd
                utc_last = utc

                if ut_start > utstart + args.offset*(utend-utstart):
                    # repeat target name if the start is delayed to make it easier to line up
                    # names and tracks
                    kwargs = {'ha' : 'right', 'va' : 'center', 'size' : int(10*args.csize)}
                    plt.text(ut_start-0.2, y, key, **kwargs)

        # Now some telescope-specific stuff (zenith holes, aerial at TNT)
        if args.telescope == 'TNT':
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

            # Plot the bad periods
            for air_start, air_end in aerial:
                plt.fill([air_start,air_end,air_end,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01], color=cols[9])
                plt.plot([air_start,air_end,air_end,air_start,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01,y-0.01], color='k',lw=0.8)

            # then for close to zenith
            start = True
            for alt, utc in zip(alts[ok], utcs[ok]):
                if start and alt > 88.:
                    air_start = utc
                    start     = False
                elif not start and alt < 88.:
                    air_end = utc
                    break

            if not start:
                plt.fill([air_start,air_end,air_end,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01], color=cols[9])
                plt.plot([air_start,air_end,air_end,air_start,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01,y-0.01], color='k',lw=0.8)

        elif args.telescope == 'VLT':
            # VLT specific
            start = True
            for alt, utc in zip(alts[ok], utcs[ok]):
                if start and alt > 86.:
                    air_start = utc
                    start = False
                elif not start and alt < 86.:
                    air_end = utc
                    break

            if not start:
                plt.fill([air_start,air_end,air_end,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01], color=cols[9])
                plt.plot([air_start,air_end,air_end,air_start,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01,y-0.01], color='k',lw=0.8)

        elif args.telescope == 'WHT' or args.telescope == 'NTT':
            # WHT / NTT specific
            start = True
            for alt, utc in zip(alts[ok], utcs[ok]):
                if start and alt > 88.:
                    air_start = utc
                    start     = False
                elif not start and alt < 87.:
                    air_end = utc
                    break

            if not start:
                plt.fill([air_start,air_end,air_end,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01], color=cols[9])
                plt.plot([air_start,air_end,air_end,air_start,air_start],
                         [y-0.01,y-0.01,y+0.01,y+0.01,y-0.01], color='k',lw=0.8)


        if key in prinfo:
            # handles the time ranges only
            pr = prinfo[key]

            pranges = pr.prange
            for t1, t2, col, lw, p_or_t in pranges:
                if p_or_t == 'Time':
                    utc1, utc2 = 24.*(t1-isun), 24.*(t2-isun)
                    if utc1 < utc_end and utc2 > utc_start:
                        ut1  = max(utc1, utc_start)
                        ut2  = min(utc2, utc_end)
                        if ut1 < ut2:
                            plt.plot([ut1,ut2],[y,y],color=cols[col],lw=lw)

        if star.eph:
            # Now the phase info

            eph = star.eph
            times = time.Time((mjd_first,mjd_last),location=site)
            if eph.time.startswith('H'):
                times += times.light_travel_time(star.position, 'heliocentric')
            elif eph.time.startswith('B'):
                times += times.light_travel_time(star.position)
            else:
                raise Exception('Unrecognised type of time = ' + eph.time)

            if eph.time == 'HJD' or eph.time == 'BJD':
                pstart = eph.phase(times[0].jd)
                pend   = eph.phase(times[1].jd)
            elif eph.time == 'HMJD' or eph.time == 'BMJD':
                pstart = eph.phase(times[0].mjd)
                pend   = eph.phase(times[1].mjd)
            else:
                raise Exception('Unrecognised type of time = ' + eph.time)

            # Compute uncertainty in predictions
            delta = 24.*min(eph.etime((pstart+pend)/2.), eph.coeff[1]/2.)
            if delta > 0.03:
                plt.plot([utend-2.*delta,utend],[y,y],'k',lw=2)
                plt.plot([utend-2.*delta,utend-2*delta],[y-0.005,y+0.005],'k',lw=2)
                plt.plot([utend,utend],[y-0.005,y+0.005],'k',lw=2)

            if key in prinfo:
                pr   = prinfo[key]

                # Draw phase ranges of interest
                pranges = pr.prange
                for p1, p2, col, lw, p_or_t in pranges:
                    if p_or_t == 'Phase':
                        d1 = pstart + (p1 - pstart) % 1 - 1
                        d2 = pend   + (p1 - pend) % 1
                        nphs = int(np.ceil(d2 - d1))
                        for n in range(nphs):
                            ut1 = utc_first + (utc_last-utc_first)*(d1 + n - pstart)/(pend-pstart)
                            ut2  = ut1 + (utc_last-utc_first)/(pend-pstart)*(p2-p1)
                            ut1  = max(ut1, utc_first)
                            ut2  = min(ut2, utc_last)
                            if ut1 < ut2:
                                plt.plot([ut1,ut2],[y,y],color=cols[col],lw=lw)

                # draws dots at phase zero
                d1 = np.ceil(pstart)
                d2 = np.floor(pend)
                nphs = int(np.ceil(d2 - d1))+1
                for n in range(nphs):
                    ut = utc_first + (utc_last-utc_first)*(d1 + n - pstart)/(pend-pstart)
                    plt.plot(ut,y,'ok')

            # draw vertical bar at meridian
            if ok.any():
                lst = mjds.sidereal_time(kind='mean', longitude=site.longitude)
                has = (lst - star.position.ra).value
                hamin, hamax = has[ok].min(), has[ok].max()
                if hamin < 0 and hamax > 0.:
                    has = np.abs(has[ok])
                    utc_mer = utcs[ok][has.argmin()]
                    plt.plot([utc_mer,utc_mer],[y-1.3*lbar,y+1.3*lbar],'k',lw=2)

    # finish off
    plt.xlabel('UTC')
    plt.title('Night starting {0!s} at the {1:s}, airmass < {2:3.1f}'.format(
            date, args.telescope, args.airmass))

    plt.xlim(utstart, utend)
    plt.ylim(0,1.05)

    axr.xaxis.set_major_locator(MultipleLocator(2))
    axr.xaxis.set_minor_locator(MultipleLocator(1))
    axr.get_xaxis().tick_bottom()
    axr.get_yaxis().set_visible(False)

    # Modify hours labels so that all lie in the range [1,24]
    fig.canvas.draw()
    labels = [item.get_text() for item in axr.get_xticklabels()]
    for i in range(len(labels)):
        try:
            if int(labels[i]) > 24:
                labels[i] = str(int(labels[i])-24)
        except:
            pass
    axr.set_xticklabels(labels)


    # add target names at left
#    pgsch(args.csize)
    for key in keys:
        y = ys[key]
        if key in prefixes:
            name = prefixes[key] + ' ' + key
        elif len(prefixes):
            name = mpre*' ' + key
        else:
            name = key
        axl.text(1,y,name,ha='right',va='center',size=int(12*args.csize))
    axl.set_xlim(0,1)
    axl.set_ylim(0,1.05)
    axl.set_axis_off()
#        pgslw(peinfo[key].lw)

    if args.hcopy:
        plt.savefig(args.hcopy)
    else:
        plt.show()


