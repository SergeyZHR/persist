#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#author Zheltoukhov Sergey 
import unittest
from astropy.io import fits
import numpy as np
import os
import persist.persist as persist

class persistTestCase(unittest.TestCase):
	def test_int(self):
		ffile =os.path.dirname(os.path.abspath(__file__) )
		persist.main('-d '+ffile+'/persist/test_fits/data/ -t '+ffile+'/tmp/ -p '+ffile+'/persist/test_fits/res/ -o')
		std=fits.open(ffile+'/persist/test_fits/res/test_dark/dark-10s-after-nlc-8.fts')
		res=fits.open(ffile+'/persist/test_fits/res0/test_dark/dark-10s-after-nlc-8.fts')
		self.assertEqual(np.nanmax(np.abs(std[0].data-res[0].data)), 0)
		std.close()
		res.close()

	def test_int_with_args(self):
		ffile =os.path.dirname(os.path.abspath(__file__) )
		persist.main_with_args(ifdat=ffile+'/persist/test_fits/data/ ',tmpdir = ffile+'/tmp/', resdir= ffile+'/persist/test_fits/res/',overwrite=True)
		std=fits.open(ffile+'/persist/test_fits/res/test_dark/dark-10s-after-nlc-8.fts')
		res=fits.open(ffile+'/persist/test_fits/res0/test_dark/dark-10s-after-nlc-8.fts')
		self.assertEqual(np.nanmax(np.abs(std[0].data-res[0].data)), 0)
		std.close()
		res.close()


	def test_name(self):
		name = persist.master_name('/azaza/test/dark-10s-after-1-if-1.fits')
		print(name)
		self.assertEqual(name, '/azaza/test/dark-10s-after-1')

	def test_if(self):
		iff = persist.ifN('/azaza/test/dark-10s-after-1-if-1.fits')
		self.assertEqual(iff, 1)



if __name__ == '__main__':
    unittest.main()
# std=fits.open('/home/szh/Documents/SciWork/persist/homb_task/res0/test/dark-10s-after-nlc-1.fts')[0].data
# res=fits.open('/home/szh/Documents/SciWork/persist/homb_task/res/test/dark-10s-after-nlc-1.fts')[0].data

# print(np.nanmean(std-res),np.nanmin(std-res))