"""

Classes for MarsSDS assembly and data retrieval

MarsSDS: Responsible for the download and creation of the
	file structure. By default, data is acquired from the station
	XB.ELYSE.

chunk, sdschunk and fdsnchunk: Classes for gathering data from
	a FDSN server or SDS modified to accept time intervals in
	Sols or UTCDateTime.

"""

from obspy.core import UTCDateTime
from obspy.core import read
from obspy.core import Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn import header

from utils import log
from utils import sol_span_in_utc
from utils import julday_in_utc
from utils import is_leap_year

from functools import reduce

import numpy as np

import os
import sys

# import pdb; pdb.set_trace()


CLIENT = Client("IRIS")


class MarsSDS(object):
	"""
	This class is responsible for the creation of the SDS as well as
	the download of the data in a given interval (in Sols)

	The time (in Sol units) is converted to UTCDateTime, so we can
	use SDSChunk to access the data.
	"""
	def __init__(self, SDSpath: str, verbose=False):

		if not os.path.exists(SDSpath):
			log(f"SDS path {SDSpath} does not exist")

		self.net      = "XB"
		self.sta      = "ELYSE"
		self.SDSpath  = SDSpath
		self.verbose  = verbose

	def __str__(self):
		return f"SetupMarsSDS:\n  Station Code: {self.net}.{self.sta}\n" +\
			   f"  SDS Path: {self.SDSpath}"

	def _make_dir(self, year, cha, data_type="D"):
		"""
		Creates the CHAN.TYPE directory in SDS if it does not exist
		"""
		year_str = str(year)
		dirs = [self.SDSpath, year_str, self.net, self.sta,
					f"{cha}.{data_type}"]
		
		for i in range(1, len(dirs) + 1):

			subdir = os.path.join(*dirs[:i])

			if not os.path.exists(subdir):
				log(f"DIR: creating {subdir}", level=0, verbose=self.verbose)
				os.mkdir(subdir)

	def _get_cha(self, loc, cha, t0, t1):
		"""
		Stores the location and channel code for all active stations
		in t0 <= t <= t1
		"""
		loc_cha_running = []

		try:
			inv = CLIENT.get_stations(network=self.net, station=self.sta,
								location=loc.replace("?", "*"),
								channel=cha.replace("?", "*"),
								starttime=t0, endtime=t1, level="channel")
		except header.FDSNNoDataException:
			code = f"{self.net}.{self.sta}.{loc}.{cha}"
			log(f"No inventory for {code} at {t0} - {t1}", level=1,
								verbose=self.verbose)
			return []

		for ch_info in inv[0][0].channels:
			# loc_cha_running["02.BHU"] = (UTCDateTime(ta), UTCDateTime(tb))
			loc_cha_running.append(f"{ch_info.location_code}.{ch_info.code}")

		return loc_cha_running

	def _gen_path(self, loc, cha, year, julday, data_type="D"):
		"""
		Returns the full SDS path of a file
		"""
		filename = ".".join([self.net, self.sta, loc, cha, data_type,
							f"{year:04d}", f"{julday:03d}"])
		return os.path.join(self.SDSpath, f"{year:04d}", self.net, self.sta,
							f"{cha}.{data_type}", filename)

	def download(self, location, channel, solstart, solend, replace=False):
		"""
		Download and store data

		location [str]: location code
		channel  [str]: channel code
		solstart [int]: initial download time in Sols
		solend   [int]: final download time in Sols
		replace [bool]: if True, will replace existing files with ney ones
		"""

		if solstart > solend or solstart < 0 or solend < 0:
			log(f"Wrong Sol: {solstart} -> {solend}")

		# Starttime and Endtime converted from Sol to UTC
		time0, _ = sol_span_in_utc(solstart)
		_, time1 = sol_span_in_utc(solend)

		# Fix continuity
		ST = Stream()

		for year in range(time0.year, time1.year + 1, 1):

			start = max(UTCDateTime(f"{year}-01-01"), time0)
			end   = min(UTCDateTime(f"{year}-12-31T23:59:59.999999"), time1)

			for julday in range(start.julday, end.julday + 1, 1):

				# Starttime and Endtime for the entire julian day
				julday_t0, julday_t1 = julday_in_utc(year, julday)

				# Active channels
				loc_cha_running = self._get_cha(location, channel, julday_t0,
												julday_t1)

				for lc_id in loc_cha_running:

					loc, cha = lc_id.split(".")
					filename = self._gen_path(loc, cha, year, julday)
					code = f"{self.net}.{self.sta}.{loc}.{cha}"

					# File already exists
					if os.path.exists(filename) and not replace:
						log(f"FILE: {os.path.basename(filename)} already " +\
									"exists. Download will skip this day",
									level=1, verbose=self.verbose)
						continue
						
					else:
						try:

							st = CLIENT.get_waveforms(network=self.net,
											station=self.sta, location=loc,
											channel=cha, starttime=julday_t0-10,
											endtime=julday_t1+10)
							ST.extend([st[0]]).merge(method=1)
							ST.trim(julday_t0-10, julday_t1+10)

							st = ST.slice(julday_t0, julday_t1,
								nearest_sample=False).select(location=loc,
								channel=cha)

							if len(st) == 0:
								continue

						except header.FDSNNoDataException as E:
							log(f"{E}", level=1, verbose=self.verbose)
							continue

						except Exception as E:
							log(f"{E}", level=2, verbose=self.verbose)

						log(f"DOWNLOAD: {code} Julday: {julday}",
											level=0, verbose=self.verbose)

						# SDS Creation/Update
						self._make_dir(year, cha)

						if any(np.ma.is_masked(tr.data) for tr in st):
							print("MASKED")
						# 	for i in range(len(st)):
						# 		st[i].data = st[i].data.filled()
						st.split().write(filename, format="MSEED")


class chunk(object):
	def __init__(self, N, S, L, C, verbose = False):
		self.N          = N
		self.S          = S
		self.L          = L
		self.C          = C
		self._verbose   = verbose

		self._S         = Stream()

	def id(self):
		return "%s.%s.%s.%s" % (self.N, self.S, self.L, self.C)

	def holdings(self):
		print(self._S)
		return

	@staticmethod
	def _syncronize(others, this):
		others.append(this)

		ss = reduce(max, map(lambda st: st[0].stats.starttime, others))
		ee = reduce(min, map(lambda st: st[0].stats.endtime, others))

		for st in others:
			st.trim(ss, ee, nearest_sample = True)

		mask = np.zeros(len(this[0]), dtype=bool)
		for st in others:
			if len(st) > 1: st.merge()
			if not hasattr(st[0].data, "mask"): continue
			mask += st[0].data.mask

		for st in others:
			st[0].data = np.ma.masked_where(mask , st[0].data)

		return

	def _clean(self, s, e = None):
		toremove = []

		for t in self._S:
			if t.stats.endtime < s:
				toremove.append(t)

		for t in toremove:
			self._S.remove(t)

		if e is not None:
			for tr in self._S:
				tr.trim(starttime=s, endtime=e)

		return

	def _update(self, s, e):
		raise NotImplemented("Please Implement Me!")

	def get(self, start, end, others, withgaps = False, incomplete = False,
						nosync = False):

		if isinstance(start, int):
			s, _ = sol_span_in_utc(start)
			_, e = sol_span_in_utc(end)

		else:
			s, e = start, end

		self._update(s, e)

		sliced_t = self._S.slice(s, e, keep_empty_traces = False,
						nearest_sample = True).copy()

		ntraces = len(sliced_t)
		sliced_t = sliced_t.merge(method = 1,
			interpolation_samples = 0,
			fill_value = None)

		if ntraces == 0: return False

		if not incomplete and (sliced_t[0].stats.starttime > s or \
						sliced_t[0].stats.endtime < e):
			log(f" W:> Trace is incomplete: {sliced_t[0].stats.starttime} >" +\
				f" {s} || {sliced_t[0].stats.endtime} < {e}.", level=1,
				verbose=self._verbose)
			return False

		if not withgaps and ntraces > 1:
			if self._verbose: print(" W:> Trace has %d gap%s." % (ntraces - 1,
						 "" if ntraces == 2 else "s"))
			return False

		if not nosync:
			self._syncronize(others, sliced_t)
		else:
			others.append(sliced_t)

		return True


class fdsnchunk(chunk):
	def __init__(self, server, N, S, L, C, user=None, password=None,
						verbose=False):
		super().__init__(N, S, L, C, verbose)

		self.server     = server
		self.user       = user
		self.password   = password
		self._start = None
		self._end   = None

		try:
			self._client = Client(server, user=user, password=password)

		except:
			log(f" No main FDSN client instantiated from {server} w/ " +\
				f"USER={user} and PASS={password}", level=2,
				verbose=self._verbose)

		log(f" I:> New FDSNCHUNK {N}.{S}.{L}.{C} @ {self.server}", level=0,
				verbose=self._verbose)

	def _update(self, s, e):
		download = self._make_download_dict(s, e)
		self._extend(download)
		return

	def _make_download_dict(self, s, e):
		s_dl, e_dl = s, e

		if e <= s:
			log(" Endtime <= Starttime: Invalid starttime/endtime", level=1,
					verbose=self._verbose)

		elif self._start is None and self._end is None:
			self._start, self._end = s, e

		elif self._start <= s and self._end >= e:
			log(f" Data from {self._start} to {self._end} was already" +\
					f" downloaded", level=0, verbose=self._verbose)
			self._start, self._end = s, e
			return None

		elif self._start <= s <= self._end <= e:
			s_dl = self._end
			self._start, self._end = s, e

		elif s < self._start:
			log(" starttime < self._start: Download will start at self._end" +\
					f" if endtime > self._end", level=1, verbose=self._verbose)
			if e > self._end:
				s_dl = self._end
				self._end = e
			else:
				self._end = e
				return None

		download = {"network"  : self.N,
					"station"  : self.S,
					"location" : self.L,
					"channel"  : self.C,
					"starttime": s_dl,
					"endtime"  : e_dl}

		return download
 
	def _extend(self, download_dict):

		if download_dict is not None:

			try:
				self._S.extend(self._client.get_waveforms(**download_dict))
				string = f"{self.N}.{self.S}.{self.L}.{self.C}:" +\
					f" {download_dict['starttime']} - {download_dict['endtime']}"
				log(f" W:> READ FDSN: {string}", level=0, verbose=self._verbose)

			except:
				log(" No data for {}.{}.{}.{} {} - {}".format(*list(download_dict.values())),
							level=2, verbose=self._verbose)

		self._S.merge(method = -1)
		self._S.sort()
		self._clean(self._start, self._end)

		return True


class sdschunk(chunk):
	"""
	This class walks a SDS archive by incremental time and
	handle the synchronization of GAPS between traces supplied.
	
	One class handle one stream, but, while the get() method is called
	a list of streams can be passed what allows for getting a syncronized
	list of traces and streams.
	
	The return of the get() method is True or False. True when data was
	added to others and False when no data was added to others.
	"""
	def __init__(self, path, N, S, L, C, verbose = False):
		super().__init__(N, S, L, C, verbose)
		
		if not os.path.isdir(path):
			log(" Bad path to SDS", level=2, verbose=verbose)
		
		self.path = path
		self.N = N
		self.S = S
		self.L = L
		self.C = C
		self._verbose = verbose
		
		self._S = Stream()
		self._last_time = None
		self._visited = []

		log(f" I:> New SDSCHUNK {N}.{S}.{L}.{C} @ {self.path}", level=0,
					verbose=verbose)

	def _make_path(self, d):
		path = "%s/%04d/%s/%s/%s.D/%s.%s.%s.%s.D.%04d.%03d" % (self.path,
					d.year, self.N, self.S, self.C,self.N, self.S, self.L,
					self.C, d.year, d.julday)

		if path in self._visited: return None
		self._visited.append(path)

		return path

	def _update(self, s, e):
		s2 = (UTCDateTime(s.date) - 86400 / 2) if UTCDateTime(s).hour == 0 and\
					UTCDateTime(s).minute < 30 else UTCDateTime(s.date)
		e2 = (UTCDateTime(e.date) + 2 * 86400) if UTCDateTime(e).hour == 23 and\
					UTCDateTime(e).minute > 30 else (UTCDateTime(e.date) + 1*86400)

		while s2 < e2:
			self._extend(self._make_path(s2))
			s2 += 86400

		self._clean(s)

		return self._last_time

	def _extend(self, filename):
		if filename is None or not os.path.isfile(filename): return False

		self._S.extend(read(filename))
		self._S.merge(method = -1)
		self._S.sort()

		log(f" W:> READ: {filename}", level=0, verbose=self._verbose)

		self._last_time = self._S[-1].stats.endtime

		return True


class DownloadFull:
	def __init__(self, N, S, L, C, output_dir):
		self.N = N
		self.S = S
		self.L = L
		self.C = C
		if not os.path.exists(output_dir):
			log(f"Output directory does not exist", level=2)
		else:
			self.output_dir = output_dir

	def download(self, start, end):
		"""
		start and end in Sols
		"""

		t0, _ = sol_span_in_utc(start)
		_, t1 = sol_span_in_utc(end)

		t0, _ = julday_in_utc(t0.year, t0.julday)
		_, t1 = julday_in_utc(t1.year, t1.julday)

		print(t0, t1)

		st = CLIENT.get_waveforms(network=self.N, station=self.S, location=self.L,
						channel=self.C, starttime=t0, endtime=t1)

		fname = f"{self.N}-{self.S}-{self.L}-{self.C}-{start}-{end}"
		st.write(os.path.join(self.output_dir, fname), format="MSEED")



if __name__ == "__main__":

	## Testing SetupMarsSDS
	# setup = MarsSDS(SDSpath="/home/rafael/mars/SDS/", verbose=True)
	# print(setup)
	# setup.download(location="02", channel="BH*", solstart=197, solend=300,
	# 					replace=True)

	# Testing SDSchunk
	# a = sdschunk("/home/rafael/mars/SDS/", "XB", "ELYSE", "02", "BHV")
	# streams = []
	# t0 = UTCDateTime("2019-05-22T23:00:00")
	# a.get(t0, t0 + 4*86400, streams, withgaps=True, incomplete=True)

	# print(len(streams))
	# import pickle
	# for st in streams:
	# 	print(st)
	# 	st.plot()
	# 	with open("/home/rafael/Desktop/stream", 'wb') as f:
	# 		pickle.dump(st, f)

	data = DownloadFull("XB", "ELYSE", "02", "BHU", "/home/rafael/mars/")
	data.download(100, 110)