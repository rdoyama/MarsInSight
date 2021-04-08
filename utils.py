from __future__ import print_function

import datetime as dt
import numpy as np
import matplotlib.mlab as mlab
import obspy
import sys

from obspy import UTCDateTime
from obspy.signal.util import next_pow_2

# Earth and Mars days in seconds
SECONDS_PER_EARTH_DAY = 86400.

# Sol0 start - Insight landing.
# Value is from the JPL URL below for Sol 0:
# https://naif.jpl.nasa.gov/cgi-bin/chronos_nsyt.pl?setup=nsyttime
REFERENCE_TIME_INSIGHT_LANDING = UTCDateTime("2018-11-26T05:10:50.336037Z")

# Sol-001 and Sol-002 start times to compute one Martian day in seconds.
# Cannot use landing time because Sol-000 lasted shorter.
SOL01_START_TIME = UTCDateTime("2018-11-27T05:50:25.580014Z")
SOL02_START_TIME = UTCDateTime("2018-11-28T06:30:00.823990Z")

# Compute one martian day in seconds with a microsecond correction
# to avoid falling into previous or next sol instead of current
# one at sol midnights
SECONDS_PER_MARS_DAY = SOL02_START_TIME - SOL01_START_TIME - 0.000005


def log(message, level=2, verbose=False):
	levels = ["Message:", "Warning:", "Error:"]
	if not hasattr(log, "lastlevel"):
		log.lastlevel = None
	
	if level == 0 and verbose is False: return
	if log.lastlevel != level:
		print("")
		print("{}".format(levels[level]))
		log.lastlevel = level
	
	print("  {}".format(message))
	if level == 2: sys.exit(1)


def sol_span_in_utc(sol, sol0_start_utc=REFERENCE_TIME_INSIGHT_LANDING):
	"""
	Returns start and end times in UTC for a given sol.
	"""
	utc_representation = \
		UTCDateTime(sol * SECONDS_PER_MARS_DAY) + float(sol0_start_utc)

	return utc_representation, utc_representation + SECONDS_PER_MARS_DAY


def julday_in_utc(year, julday):
	"""
	Returns start and end times in UTC for a given julian day
	"""
	julday1 = UTCDateTime(f'{year}-01-01')
	julday1_end = UTCDateTime(f'{year}-01-01T23:59:59.999999')

	julday_s = julday1 + 86400.0 * (julday - 1)
	julday_e = julday1_end + 86400.0 * (julday - 1)

	return julday_s, julday_e

def julday_in_sol(year, julday):
	"""
	Converts Julian day to InSight Mission Sol
	Returns the Sol at the start and at the end of the julian day,
	respectively.
	"""
	julday_s, julday_e = julday_in_utc(year, julday)

	if julday_s < REFERENCE_TIME_INSIGHT_LANDING:
		log(f"Day {year:04d}.{julday:03d} before the start of the mission")

	if REFERENCE_TIME_INSIGHT_LANDING <= julday_s <= SOL01_START_TIME:
		return 0, 1

	sol1 = int((julday_s - SOL01_START_TIME) / SECONDS_PER_MARS_DAY)
	sol2 = int((julday_e - SOL01_START_TIME) / SECONDS_PER_MARS_DAY)

	return sol1, sol2


def is_leap_year(year):
	"""
	Returns True if "year" is a leap year, else False
	"""
	if year % 4 != 0:
		return False
	elif year % 100 != 0:
		return True
	elif year % 400 != 0:
		return False
	else:
		return True

