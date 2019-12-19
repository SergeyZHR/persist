#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Программа работает со всеми файлами в каталоге (имя файла без пути не длиннее 120 символов!)
# выбираются все sif файлыБ разбиваются на соответствующие группы и корректируются.
# результат пишется в файл с таким же базовым именем, но с добавкой nlc
# Коррекция за нелинейность по готовым коэффициентам параболы - файлы c_coef_master.fts, b_coef_master.fts, a_coef_master.fts (коэф. при x2)
# и готовому файлу bias_master (используется вычисляемый файл) - файл не влияет на конечный ответ, куда выводится полный накопленный сигнал
# все файлы коэф. и bias готовятся программой mnk2.py

import sys
import getpass
import numpy as np
import os
import re
#import time
import shlex
import argparse
import persist.overlight as overlight
#import badproc

#мой модуль для сортировки ANC NDR-файлов 
from persist.list_files import sort_files as sort_files

#стандартный модуль для подсчета в массиве одинаковых элементов
from collections import Counter

from datetime import datetime
import glob

#from ctypes import *
from astropy.io import fits
from astropy.time import Time
#from astropy.modeling import models, fitting
#from scipy.ndimage import median_filter

def main(string):

	ndr_num_change=False

	#формирует куб данных
	def get_imagedata(Imagelist):
		global ndr_num_change
		for i, sImage in enumerate(Imagelist):
			hdulist = fits.open(sImage)
			if i==0:
				imagedata = np.array(hdulist[0].data, dtype=np.float32)
				if len(Imagelist)==1:
					imagedata = np.tile(imagedata,(1,1))
				else:
					imagedata = np.tile(imagedata,(len(Imagelist),1,1))
			else:
				imagedata[i,:,:] = np.array(hdulist[0].data, dtype=np.float32)
		itime=hdulist[0].header['ITIME']
		ndrs=hdulist[0].header['NDRS']
		if ndrs >= ndr_min:
			ndr_num_change = True
			imagedata = imagedata[1:,:,:]
		else:
			ndr_num_change = False
		xx=hdulist[0].header['X']
		yy=hdulist[0].header['Y']
		wid=hdulist[0].header['WID']
		hgt=hdulist[0].header['HGT']
		prihdr = hdulist[0].header
		hdulist.close()
		
		return (imagedata,itime,ndrs,xx,yy,wid,hgt,prihdr)

	def createParser():
		parser = argparse.ArgumentParser(prog='non_lin_corr7.py', description ='non linear correction')

		parser.add_argument ('-d', '--directory', type = str, 
			help = 'Directory with files (fits, fit, and fts) which must be corrected (only one)')

		parser.add_argument ('-p', '--distanationPoint', type = str, 
			help = 'Directory for result frames. If not specified, will be /home/anc/data/results/')

		parser.add_argument ('-o', '--overwrite', action='store_const', const=True, default=False, 
			help = 'Overwrite already reduced files. If not specified, programm will skip them.')

		parser.add_argument ('-l', '--overlight', action='store_const', const=True, default=False, 
			help = 'correct overlight pixels')

		parser.add_argument ('-S', '--silent', action='store_const', const=True, default=False,
			help = 'prevent do not display a lot of information')

		parser.add_argument ('-N', '--minNDRS', type=int, 
			help = 'the minimum number of NDRs at which the first is not taken into account')

		parser.add_argument ('-b', '--bad', action='store_const', const=True, default=False, 
			help = 'correct bad pixels')



		return parser

	parser = createParser()

	namespace = parser.parse_args(shlex.split(string))

	#путь к коэф. NLC
	patch=os.path.dirname(os.path.abspath(__file__) )+'/coef/'

	#выбираем куда пишем результат
	if namespace.distanationPoint == None:
		end_patch = '/home/anc/data/results/'
	else:
		if namespace.distanationPoint[-1] != '/' : namespace.distanationPoint+='/'
		end_patch = namespace.distanationPoint


	if namespace.directory[-1] != '/' : namespace.directory+='/'
	file_patch = namespace.directory

	if not os.path.exists(end_patch):
		os.makedirs(end_patch)
		pass

	is_overlight_correct = namespace.overlight   #OVERLIGHT
	is_bad_correct = namespace.bad   #bad pixels

	if namespace.minNDRS!= None:
		ndr_min = namespace.minNDRS
		if ndr_min <3:
			if not namespace.silent:print('algorithm need at least 2 ndrs to work, minNDRS  = 3')
			ndr_min = 3
	else:
		ndr_min = np.inf

	#
	#txt_hot_exists=0
	#try:
	#	y_hot,x_hot=np.loadtxt(patch+'hot_pixels.dat', dtype='int16', usecols=[0,1], unpack=True)
	#except IOError:
	#	print 'Hot pixels data file is absent'
	#else:
	#	txt_hot_exists=1
	#	print "Hot pixels data file exists"


	vstavka='#*#'
	a_coef_med=4.533e-12
	b_coef_med=2.517e-6
	c_coef_med=-0.00686
	bias_med=58101.14

	# список вообще всех sif файлов
	list_all_sif=glob.glob(file_patch+'*-if-*.fits')

	# вспомогательный массив имен
	bna=np.zeros(len(list_all_sif),dtype='|S120')

	for i in range(len(list_all_sif)):
		sp=list_all_sif[i][0:list_all_sif[i].rfind("-if")]
		bna[i]=sp[0:sp.rfind("-")]

	# список всех базовых имен
	base_name_all=list(Counter(bna))
	#цикл перебора всех базовых имен
	for kkk in range(len(base_name_all)):
		base_name=base_name_all[kkk].decode()
		if not namespace.silent:print(base_name)

		#l=glob.glob(base_name+"-?-if*.fits")+glob.glob(base_name+"-??-if*.fits")+glob.glob(base_name+"-???-if*.fits")+glob.glob(base_name+"-????-if*.fits")
		l=glob.glob(base_name+"-*-if*.fits")

		num=np.zeros(len(l),dtype=np.int)

		for i in range(len(l)):
			sp=l[i][0:l[i].rfind("-if")]
			num[i]=int(sp[sp.rfind("-")+1:])

	#список уникальных номеров ramp'ов с базовым именем base_name
		count=list(Counter(num))
	#делаем цикл по списку count
	# для каждого номера формируем новое базовое имя - собираем все NDR для него
	# делаем коррекцию за нелинейность и считаем ramp

		for iii in range(len(count)):
			try:
				if not namespace.silent:print("	START",datetime.now() )
				if not namespace.silent:print("		FileName", base_name+'-'+str(count[iii])+'-if*.fits')
				if not namespace.silent:print("		len(count)", len(count))
				list_of_images=sort_files(glob.glob(base_name+'-'+str(count[iii])+'-if*.fits'))

		#reading TXT file with OCS data and preparing KEYS-VALUES for FITS header
				txt_exists=0
				try:
					txt_file=open(base_name+'-'+str(count[iii])+'.txt')
				except IOError:
					if not namespace.silent:print('OCS data (TXT '+base_name+'-'+str(count[iii])+'.txt'+') file opening error')
				else:
					txt_exists=1
					if not namespace.silent:print("TXT file exists")
					ocs_data=txt_file.readline()
					ocs_data.rstrip('\n')
					#replace all {} with NONE
					ocs_data=re.sub(r'{}', 'NONE', ocs_data)
					#search all values with space and replace it with vstavka
					list_bad_values=re.findall(r'\{[^\{]{,100}\}', ocs_data)
					if not namespace.silent:print("list", list_bad_values)
					for bv in list_bad_values:
						bvv=bv.replace(' ',vstavka)
						ocs_data=ocs_data.replace(bv, bvv)
					d=dict()
					ocs_data=ocs_data.replace('{','')
					ocs_data=ocs_data.replace('}','')
					ocs_data2=ocs_data.split(' ')
					for i in range(0, len(ocs_data2), 2):
						d[ocs_data2[i]]=ocs_data2[i+1]
				
				#data - куб данных
				(data, ITIME, NDRS, XX, YY, WID, HGT, HDR) = get_imagedata(list_of_images)

				if not namespace.silent:print("ITIME=", ITIME, "NDRS=", NDRS, "X=", XX, "Y=", YY, "WID=", WID, "HGT=", HGT)
		# dTIME - время между 2-мя последовательными NDR

				dTIME=ITIME/(NDRS-1)

				x = np.array([i for i in range(data.shape[0])])
				x=(x+1)*dTIME

				if not namespace.silent:print("3D array ready",datetime.now())

		# proverka togo, 4to kadr lezhit v rabo4ej oblasti.
		# esli eto tak, to ispol'zuem individ. coef. non lin. corr.
		# esli net, to dlja vsego kadra berem mediannye

				if (XX>=512) and (YY>=512) and (XX+WID<=1536) and (YY+HGT<=1536) :

		#читаем bias и коэф. корректирующей параболы (размер области, где они определены 512,512+1024x1024)
					
					if not namespace.silent:print("Individual correction")
					
					hdulist=fits.open(patch+'bias_master.fts')
					bias_full=hdulist[0].data

					hdulist=fits.open(patch+'a_coef_master.fts')
					a_coef_full=hdulist[0].data

					hdulist=fits.open(patch+'b_coef_master.fts')
					b_coef_full=hdulist[0].data

					hdulist=fits.open(patch+'c_coef_master.fts')
					c_coef_full=hdulist[0].data

		#использовать часто надо не всю область, т.к. снимки были получены с части приемника
		# выделяем и присваиваем нужную часть
					a_coef=a_coef_full[YY-512:YY-512+HGT, XX-512:XX-512+WID]
					b_coef=b_coef_full[YY-512:YY-512+HGT, XX-512:XX-512+WID]
					c_coef=c_coef_full[YY-512:YY-512+HGT, XX-512:XX-512+WID]
					bias=bias_full[YY-512:YY-512+HGT, XX-512:XX-512+WID]
				elif (XX<=512) and (YY==0) and (WID>=1024) and (HGT==2048) :
					hdulist=fits.open(patch+'bias_master.fts')
					bias_full=hdulist[0].data

					hdulist=fits.open(patch+'a_coef_master.fts')
					a_coef_full=hdulist[0].data

					hdulist=fits.open(patch+'b_coef_master.fts')
					b_coef_full=hdulist[0].data

					hdulist=fits.open(patch+'c_coef_master.fts')
					c_coef_full=hdulist[0].data
					a_coef=np.zeros((HGT,WID), dtype=np.float32)+a_coef_med
					b_coef=np.zeros((HGT,WID), dtype=np.float32)+b_coef_med
					c_coef=np.zeros((HGT,WID), dtype=np.float32)+c_coef_med
					bias=np.zeros((HGT,WID), dtype=np.float32)+bias_med
					#a_coef[512:512+1024,512-XX:512-XX+1024]=a_coef_full
					#b_coef[512:512+1024,512-XX:512-XX+1024]=b_coef_full
					#c_coef[512:512+1024,512-XX:512-XX+1024]=c_coef_full
					#bias[512:512+1024,512-XX:512-XX+1024]=bias_full
					
					#отступили 2 пиксела от края
					a_coef[517:512+1019,517-XX:512-XX+1019]=a_coef_full[5:1019,5:1019]
					b_coef[517:512+1019,517-XX:512-XX+1019]=b_coef_full[5:1019,5:1019]
					c_coef[517:512+1019,517-XX:512-XX+1019]=c_coef_full[5:1019,5:1019]
					bias[517:512+1019,517-XX:512-XX+1019]=bias_full[5:1019,5:1019]
				else:
					if not namespace.silent:print("Average (median) correction")
					a_coef=a_coef_med
					b_coef=b_coef_med
					c_coef=c_coef_med
					bias=bias_med
					
					#warm pixel processing 
				if is_bad_correct: 
					bad_x,bad_y,bad_value = badproc.processing_cool_warm_pixels(data,XX,YY,WID,HGT,ITIME)

				for m in range(data.shape[0]):
					y=bias-data[m,:,:]
					data[m,:,:]=bias-y*(1+c_coef+b_coef*y+a_coef*y**2)


				if not namespace.silent:print("Begin", datetime.now())
				n=data.shape[0]

				sum_x=0.0
				sum_x2=0.0

				for i in range(x.shape[0]):
					sum_x=sum_x+x[i]
					sum_x2=sum_x2+x[i]**2

				sum_y=np.zeros((data.shape[1],data.shape[2]),dtype=np.float32)
				sum_xy=np.zeros((data.shape[1],data.shape[2]),dtype=np.float32)

				for m in range(data.shape[0]):
					sum_y=sum_y+data[m,:,:]
					sum_xy=sum_xy+data[m,:,:]*x[m]

				znam=n*sum_x2-sum_x**2

				coef_a=np.zeros((data.shape[1],data.shape[2]),dtype=np.float32)
				coef_b=np.zeros((data.shape[1],data.shape[2]),dtype=np.float32)

				coef_b=(sum_y*sum_x2-sum_x*sum_xy)/znam
				coef_a=-ITIME*(n*sum_xy-sum_x*sum_y)/znam

				med_frame=np.nanmedian(coef_a.flatten())
				coef_a[np.where(coef_a > 9.e4)]=med_frame


				#WARM PIXELS
				#coeff_a[] = warm.warm()

				if is_overlight_correct:	   #OVERLIHT
					tt=bias
					zeroLevel = bias-tt*(1+c_coef+b_coef*tt+a_coef*tt**2)
					coef_a = overlight.correct(data=data, res=coef_a, points=x, zeroLevel=zeroLevel, scale=-ITIME, bias = bias,nnc=ndr_num_change )	 #OVERLIHT
			

		#hot pixels correction
				if is_bad_correct:
					if not namespace.silent:print("bad pixel correction start")
					coef_a[bad_y,bad_x] = bad_value
					coef_a = badproc.processing_bad_pixels(coef_a,XX,YY,WID,HGT)

					if not namespace.silent:print("bad pixel correction end")
				

		#coef_a - сигнал
				hdu=fits.PrimaryHDU(coef_a)
				hdu.header=HDR

		# формируем DATE-OBS в FITS Header
				s_date=hdu.header['DATE']
				ss=s_date.split('.')
				a_date='-'
				aa=a_date.join(ss)+'T'+hdu.header['TIME']
				times=Time(aa)
				hdu.header.set('DATE-OBS',times.value)
				hdu.header.set('BZERO',0)

				if hdu.header['UPPER'] !='OPEN' and hdu.header['LOWER']!= 'OPEN':
					hdu.header.insert(45,('FILTER',hdu.header['UPPER']+'_'+hdu.header['LOWER']))
				else:
					if hdu.header['UPPER'] != 'OPEN' :
						hdu.header.insert(45,('FILTER',hdu.header['UPPER']))
					if hdu.header['LOWER'] != 'OPEN' :
						hdu.header.insert(45,('FILTER',hdu.header['LOWER']))
						
				hdu.header.set('EXPTIME',hdu.header['ITIME'])

				if txt_exists==1:
		# place OCS data in FITS header
					if d['CURDEC'].find('-') > -1:
						d['CURDEC']='-'+d['CURDEC'].replace('-','')

					if d['TARDEC'].find('-') > -1:
						d['TARDEC']='-'+d['TARDEC'].replace('-','')

					for s_tmp in d.keys():
						if len(s_tmp)<9:
							try:
								dig1=int(d[s_tmp].rstrip())
								hdu.header.set(s_tmp,dig1)
							except:
								try:
									dig2=float(d[s_tmp].rstrip())
									hdu.header.set(s_tmp,dig2)
								except:
									hdu.header.set(s_tmp, d[s_tmp].rstrip().replace(vstavka, ' '))
							

				try:
					hdu.writeto(end_patch+base_name.split('/')[-1]+'-nlc-'+str(count[iii])+'.fts',overwrite=namespace.overwrite)
				except:
					print('error write result!')

				if not namespace.silent:
					print("	END", datetime.now())
				else:
					print(base_name.split('/')[-1]+'-nlc-'+str(count[iii])+'.fts NLC corrected')

			except KeyboardInterrupt:
				sys.exit(1)

		# 	except Exception as e:
		# 		print('\n\n\n')
		# 		print('nlc_ERROR!  '+ file_patch)
		# 		print(str(e))
		# 		print('\n\n\n')
	 # 