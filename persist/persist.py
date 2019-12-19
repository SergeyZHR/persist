#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#author Zheltoukhov Sergey 
from astropy.io import fits
import numpy as np
import glob	
import datetime
import pandas as pd
from subprocess import call
import os
import sys
import argparse
import time
import persist.my_non_lin_corr7 as my_non_lin_corr7

sss=time.time()

a0 = 	[0.2,	0.7,	0.7,	0.9,	]
alph0 = [14,	3.5,	0.55,	0.09,	]
b0 = 	[2.15,	1.9,	1.5,	1.7		]
bet0 = 	[-1,	0.2,	0.11,	0.025	]
a1max = 35
t01=17
t02 = 90
t03 = 500
t04 = 4800


def createParser():
	parser = argparse.ArgumentParser(prog='persist_correct.py', description ='persistance correction')

	parser.add_argument ('-d', '--ifdat', type = str, 
		help = 'Date directory with raw, -if- files ')

	parser.add_argument ('-n', '--nlcdat', type = str, 
		help = 'Date directory with nlc files, if not specified program replace "data" to "data_nlc" in ifdat ')

	parser.add_argument ('-t', '--tmpdir', type = str, 
		help = 'Directory for temporary data, if not specified, program use /data/persist/tmp/ ')

	parser.add_argument ('-p', '--resdir', type = str, 
		help = 'Directory for result data, if not specified, program try to overwrite nlc files! ')

	parser.add_argument ('-s', '--subfolders', type = str, nargs='*',
		help = 'subfolders in -if- directory for processing, if not specified, program process all sub(sub...)folders ')

	parser.add_argument ('-o', '--overwrite', action='store_const', const=True, default=False, 
		help = 'Overwrite already reduced files. If not specified, programm will skip them.')



	return parser


def fa(x,a,alpha,b,beta):
	#x1=np.array(x)      # asarray 
	x1=np.asarray(x)
	ans = alpha*(x1-a)	
	x2 = x1[x1>b]		
	ans[x1>b] = beta*(x2-b)+alpha*(b-a)
	return ans

def persistA(l_raw):
	l = l_raw/60000

	l[l>4]=4

	a1 = fa(l,a0[0],alph0[0],b0[0],bet0[0])
	a1[a1<0]=0
	a1[a1>a1max]=a1max

	a2 = fa(l,a0[1],alph0[1],b0[1],bet0[1])
	a2[a2<0]=0

	a3 = fa(l,a0[2],alph0[2],b0[2],bet0[2])
	a3[a3<0]=0

	a4 = fa(l,a0[3],alph0[3],b0[3],bet0[3]) 
	a4[a4<0]=0

	return np.dstack((a1,a2,a3,a4))

def persistT(dt):
	return np.array([np.exp(-1*dt/t01), np.exp(-1*dt/t02), np.exp(-1*dt/t03), np.exp(-1*dt/t04) ])#.astype(np.float)  #dtype float


def ifN(s):
	if 'if' not in s:
		return 0
	else:
		return int(s[s.find('-',-9)+1:][:-5])

def master_name(s):
	if 'if' not in s:
		return -1
	else:
		return s[:s.find('-',-9)+1][:-4]



def main(string=None):

	parser = createParser()
	if string!=None:
		import shlex
		namespace = parser.parse_args(shlex.split(string))
	else:	
		namespace = parser.parse_args()


	path=os.path.dirname(os.path.abspath(__file__) )+'/'

	# nlc_dat = '/home/sergey/Documents/SciWork/persist/data_nlc/20190722/'
	# if_dat = '/home/sergey/Documents/SciWork/persist/data/20190722/'

	# temp_dir = '/home/sergey/Documents/SciWork/persist/tmp/'
	# res_dir = '/home/sergey/Documents/SciWork/persist/res/'


	if namespace.tmpdir == None:
		temp_dir = '/data/persist/tmp/'
	else:
		if namespace.tmpdir[-1] != '/' : namespace.tmpdir+='/'
		temp_dir = namespace.tmpdir

	if namespace.ifdat[-1] != '/' : namespace.ifdat+='/'
	if_dat = namespace.ifdat

	if namespace.nlcdat == None:
		nlc_dat = if_dat.replace('data','data_nlc')
	else:
		if namespace.nlcdat[-1] != '/' : namespace.nlcdat+='/'
		nlc_dat = namespace.nlcdat

	if namespace.resdir == None:
		res_dir = nlc_dat
	else:
		if namespace.resdir[-1] != '/' : namespace.resdir+='/'
		res_dir = namespace.resdir


	listFF=[]
	times = []
	exp = []
	ndrs = []
	xx = []
	yy = []
	wid = []
	hgt = []

	dirs_nlc  = []
	for path_t, _, _ in os.walk(nlc_dat):
		if path_t != nlc_dat:
			dirs_nlc.append(path_t+'/')
		else:
			dirs_nlc.append(path_t)

	for directory in dirs_nlc:
		listNLC = glob.glob(directory+'*nlc*.fts')
		for file in listNLC:
			fitsFile = fits.open(file)
			listFF.append(file)
			times.append(datetime.datetime.strptime(fitsFile[0].header['DATE']+' '+fitsFile[0].header['TIME'],'%Y.%m.%d %H:%M:%S.%f') )
			exp.append(fitsFile[0].header['ITIME'] )
			ndrs.append(fitsFile[0].header['NDRS'] )
			xx.append(fitsFile[0].header['X'] )
			yy.append(fitsFile[0].header['Y'] )
			wid.append(fitsFile[0].header['WID'] )
			hgt.append(fitsFile[0].header['HGT'] )

			fitsFile.close()

	dataILUM = pd.DataFrame({'file': listFF, 'time_0': times, 'exp':exp , 'ndrs':ndrs, 'X':xx, 'Y':yy, 'WID':wid, 'HGT':hgt})
	dataILUM['time'] = dataILUM.apply(lambda t : t['time_0'] + datetime.timedelta(seconds = t['exp'] ) ,axis = 1)
	dataILUM=dataILUM.sort_values('time').reset_index()
	dataILUM = dataILUM.drop('index',1)
	dataILUM['name'] = dataILUM['file'].map(lambda s: s.split('/')[-1])
	dataILUM['is_add'] = False


	listFF=[]
	times = []
	exp = []
	ndrs = []
	xx = []
	yy = []
	wid = []
	hgt = []

	dirs_if  = []
	if namespace.subfolders == None:
		for path_t, _, _ in os.walk(if_dat):
			if path_t != if_dat:
				dirs_if.append(path_t+'/')
			else:
				dirs_if.append(path_t)
	else:
		namedirs = glob.glob(if_dat+'*/')
		for namedir in namedirs:
			for subfolder in namespace.subfolders:
				if namedir.split('/')[-2] == subfolder:
					dirs_if.append(namedir)
	if len(dirs_if) ==0:
		print('no folders to correct')
		exit()

	for directory in dirs_if:
		listIF = glob.glob(directory+'*-if-*.fits')
		for file in listIF:
			fitsFile = fits.open(file)
			listFF.append(file)
			times.append(datetime.datetime.strptime(fitsFile[0].header['DATE']+' '+fitsFile[0].header['TIME'],'%Y.%m.%d %H:%M:%S.%f') )
			exp.append(fitsFile[0].header['ITIME'] )
			ndrs.append(fitsFile[0].header['NDRS'] )
			xx.append(fitsFile[0].header['X'] )
			yy.append(fitsFile[0].header['Y'] )
			wid.append(fitsFile[0].header['WID'] )
			hgt.append(fitsFile[0].header['HGT'] )
			fitsFile.close()
		
	dataIF = pd.DataFrame({'file': listFF, 'time_0': times, 'exp':exp , 'ndrs':ndrs, 'X':xx, 'Y':yy, 'WID':wid, 'HGT':hgt})
	dataIF['ifN'] = dataIF['file'].map(ifN)
	dataIF['time'] = dataIF.apply(lambda t: t['time_0'] + datetime.timedelta(seconds = t['exp'] * t['ifN']/t['ndrs'] ), axis = 1 )
	dataIF['name'] = dataIF['file'].map(lambda s: s.split('/')[-1])
	dataIF=dataIF.sort_values('time').reset_index()
	dataIF = dataIF.drop('index',1)

	targets = dataIF[dataIF['ifN']>0]#.loc[170:210]#.loc[:'dark-10s-after-5-if-2.fits']#.iloc[::2]#.loc['dark-10s-after-nlc-4.fts':'dark-10s-after-nlc-5.fts'] 
	targets['master_name'] = targets['file'].map(master_name)

	#print(targets)

	targetGroups = targets.groupby('master_name', sort = False)

	a_arr = persistA(np.zeros((2048,2048) ) ).astype(np.float)
	last_time = dataILUM.loc[0,'time']

	folders_corr = []
	tmp_ifs = []

	for tar_name, tar_ifs in targetGroups:
	#	print(tar_name)
		tar_ifs=tar_ifs.reset_index()
	#	print(tar_ifs['ifN'])

		for index, frame in dataILUM.iterrows():
			if tar_ifs.loc[0,'time_0'] <= frame['time_0']+datetime.timedelta(seconds = 0.5 ) :
				continue
			try:
				fitsFile = fits.open(frame['file'])
				l = np.zeros((2048,2048))
				# print(type(frame['X']))
				# print(frame['X'],frame['Y'],frame['HGT'],frame['HGT'])
				l[frame['Y']:frame['Y']+frame['HGT'],frame['X']:frame['X']+frame['WID'] ] =  fitsFile[0].data
				d_t0 = (frame['time']-last_time).total_seconds()
				da_arr = persistA(l)
				a_arr = a_arr * persistT(d_t0) + da_arr 
				last_time = frame['time']
				dataILUM.loc[index,'is_add']=True
				print('persistence from',frame['name'], 'add to array')

			except IOError as e:				#fix err
				print(str(e))
				print('error on',frame['name'])
				dataILUM.loc[index,'is_add']=True
				pass 
			sys.stdout.flush()
		dataILUM = dataILUM[~dataILUM['is_add']].reset_index().drop('index',1)			


		pers = np.zeros((2048,2048)).astype(np.float)
		for i in range(len(tar_ifs)):
			try:
				target=tar_ifs.iloc[i]

				dt_if = (target['time'] - last_time).total_seconds()

				pers += np.sum(a_arr * persistT(dt_if), axis = 2 ) * target['exp'] /(target['ndrs'] - 1. )
				# print(np.nanmax(pers[target['Y']:target['Y']+target['HGT'],target['X']:target['X']+target['WID'] ]))

				tar_fits = fits.open(target['file'])
				tar_fits[0].data = tar_fits[0].data.astype(float) + pers[target['Y']:target['Y']+target['HGT'],target['X']:target['X']+target['WID'] ]

				res_name = temp_dir + target['file'][len(if_dat):]
				os.makedirs(os.path.dirname(res_name),exist_ok=True)
				folders_corr.append(os.path.dirname(res_name))
				tmp_ifs.append(res_name)
				tar_fits.writeto(res_name,overwrite=namespace.overwrite)
				tar_fits.close()
				print(target['name'], 'persist corrected')

			except Exception as e:
				print(str(e))
				print('Error on', target['name'])
			sys.stdout.flush()


	print('time,', time.time() -sss )

	folders_corr = np.unique(np.array(folders_corr))

	for folder in folders_corr:
		res_folder = res_dir + folder[len(temp_dir):]
		if namespace.overwrite:
			over = ' -o'
		else:
			over = ''
#		call('/home/szh/Documents/SciWork/spectra/nlc_pipe/my_non_lin_corr.py -lb -d '+folder+'  -S -p '+res_folder+over ,shell = True)
		my_non_lin_corr7.main( '-l -d '+folder+'  -S -p '+res_folder+over )


	for tfile in tmp_ifs:
		try:
			os.remove(tfile)			
		except:
			pass

	for folder in folders_corr:
		try:
			os.removedirs(folder)
		except:
			pass

if __name__=='__main__':
	main()