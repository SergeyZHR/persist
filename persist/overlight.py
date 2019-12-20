#!/usr/bin/python3
# -*- coding: utf-8 -*-
#from astropy.io import fits
#author Zheltoukhov Sergey 
import numpy as np
#import time
#import os
#import glob


def correct(data,res,points,zeroLevel,scale,bias,nnc=False):
	"""поиск пересвеченных, ушедших в 0, в процессе наблюдений, и их исправление"""

	lastF = data[-1]

	#print points
	dADU = 100
	overPixels = np.where(lastF<zeroLevel+dADU)  	#ищем пиксели, нулевые в последнем считывании
	
	count = 0
	
	for i in range(len(overPixels[0])):			#цикл по найденным пикселям

		x,y = [ overPixels[0][i],overPixels[1][i] ]
		pixel = data[:,x,y]						#все считывания этого пикселя
		try:
			last_norm = ( np.nonzero(pixel>zeroLevel[x,y]+dADU) )[0][-1]
		except IndexError as e:
			#print(0)
			last_norm = -1

		if last_norm>0:				#хотим что бы было хотя бы 2 не нулевых считывания
			count+=1
			norm_pixel = data[:last_norm+1,x,y]
			#print data[:,x,y],norm_pixel, zeroLevel[x,y]
			value = np.polyfit(points[:last_norm+1],norm_pixel ,1)[0]  	#считаем МНК по ненулевым считываниям
			# print(norm_pixel)
			# print(data[:,x,y])
			# print(res[x,y], np.polyfit(points, data[:,x,y] ,1)[0]*scale, value*scale)
			res[x,y] = value*scale			#заменяем элемент

		if last_norm==0:				#если есть только одно
			count+=1
			value = (data[0,x,y]-bias[x,y])/points[0]  	#считаем отклонение от bias
			# print(norm_pixel)
			# print(data[:,x,y])
			#print(res[x,y], value*scale,bias[x,y],data[0,x,y],data[1,x,y],zeroLevel[x,y])
			res[x,y] = value*scale			#заменяем элемент
		if last_norm==-1:			#если вообще нет
			if nnc:
				res[x,y] = 65536*(len(points)-1)
			else:						
				res[x,y] = 65536*(len(points))/2.


	print (count, 'overlight pixel corrected')
	return res

		#print (last_norm)
