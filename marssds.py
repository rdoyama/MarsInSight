# Class for MarsSDS assembly and data retrieval

from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

from utils import log
from utils import sol_span_in_utc
from utils import julday_in_utc
from utils import is_leap_year

import os

# import pdb; pdb.set_trace()


CLIENT = Client("IRIS")


class SetupMarsSDS:
	"""
	This class is responsible for the creation of the SDS as well as
	the download of the data in a given interval (in Sols)

	The time (in Sol units) is converted to UTCDateTime, so we can
	use SDSChunk to access the data.
	"""
	def __init__(self, loc: str, cha: str, SDSpath: str, solstart: int,
					solend: int, verbose=False):

		if not os.path.exists(SDSpath):
			log(f"SDS path {SDSpath} does not exist")

		if solstart > solend or solstart < 0 or solend < 0:
			log(f"Wrong Sol: {solstart} -> {solend}")

		self.net      = 'XB'
		self.sta      = 'ELYSE'
		self.loc      = loc
		self.cha      = cha
		self.SDSpath  = SDSpath
		self.solstart = solstart
		self.solend   = solend
		self.verbose  = verbose

	def __str__(self):
		code = f"{self.net}.{self.sta}.{self.loc}.{self.cha}"
		t0, _ = sol_span_in_utc(self.solstart)
		_, t1 = sol_span_in_utc(self.solend)

		return f"SetupMarsSDS:\n  Code: {code}\n  Start: Sol {self.solstart}" +\
			   f" = {t0}\n  End: Sol {self.solend} = {t1}"

	def _make_dir(self, year: int, cha: str, data_type='D'):
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

	def _get_cha(self, t0, t1):
		"""
		Stores the location and channel code for all active stations
		in t0 <= t <= t1
		"""
		loc_cha_running = []

		inv = CLIENT.get_stations(network=self.net, station=self.sta,
							location=self.loc.replace("?", "*"),
							channel=self.cha.replace("?", "*"),
							starttime=t0, endtime=t1, level='channel')

		for ch_info in inv[0][0].channels:
			# loc_cha_running["02.BHU"] = (UTCDateTime(ta), UTCDateTime(tb))
			loc_cha_running.append(f"{ch_info.location_code}.{ch_info.code}")

		# log(f"A total of {len(loc_cha_running)} channels are working",
		# 					level=0, verbose=self.verbose)

		return loc_cha_running

	def _gen_path(self, loc, cha, year, julday, data_type='D'):
		"""
		Returns the full SDS path of a file
		"""
		filename = '.'.join([self.net, self.sta, loc, cha, data_type,
							f"{year:04d}", f"{julday:03d}"])
		return os.path.join(self.SDSpath, f"{year:04d}", self.net, self.sta,
							f"{cha}.{data_type}", filename)

	def download(self, replace=False):
		"""
		Download and store data

		replace: if it is True, will replace existing files with ney ones
		"""

		# Starttime and Endtime converted from Sol to UTC
		time0, _ = sol_span_in_utc(self.solstart)
		_, time1 = sol_span_in_utc(self.solend)


		# print(f"Years: {time0.year} - {time1.year}")

		for year in range(time0.year, time1.year + 1, 1):

			# print(year)

			start = max(UTCDateTime(f"{year}-01-01"), time0)
			end   = min(UTCDateTime(f"{year}-12-31T23:59:59.999999"), time1)

			# print(start, end)

			for julday in range(start.julday, end.julday + 1, 1):
				julday_t0, julday_t1 = julday_in_utc(year, julday)

				# Active channels
				loc_cha_running = self._get_cha(julday_t0, julday_t1)

				for lc_id in loc_cha_running:
					loc, cha = lc_id.split('.')

					# Starttime and Endtime for the entire julian day
					# print("###")
					# print(julday_t0, cha_t0)
					# print(julday_t1, cha_t1)
					# print("###")
					# download_s = max(julday_t0, cha_t0)
					# download_e = min(julday_t1, cha_t1)

					# print(download_s, download_e)

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
											channel=cha, starttime=julday_t0,
											endtime=julday_t1)

						except FDSNNoDataException as E:
							log(f"{E}", level=1, verbose=self.verbose)
							continue

						except Exception as E:
							log(f"{E}", level=2, verbose=self.verbose)

						log(f"DOWNLOAD: {code} Julday: {julday}",
											level=0, verbose=self.verbose)

						# SDS Creation/Update
						self._make_dir(year, cha)

						st.write(filename, format='MSEED')


if __name__ == '__main__':
	setup = SetupMarsSDS(loc='02', cha='BH*', SDSpath='/home/rafael/mars/SDS/',
						solstart=388, solend=391, verbose=True)
	print(setup)
	setup.download(replace=True)

