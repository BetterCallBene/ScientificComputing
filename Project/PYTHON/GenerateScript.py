#!/usr/bin/python

import os
import re
import numpy as np
from shutil import move


#MAPLE_PATH = '/Library/Frameworks/Maple.framework/Versions/Current/bin/maple'
#CONFIG_NAME = 'config.py'

def generate(MAPLE_PATH, template_path, output, maple_path, tmp_files, reps):
	
	print 'Temporaere loeschen'
	
	for filename in tmp_files:
		if os.path.exists(filename):
			os.remove(filename)

	if os.path.exists(output):
		os.remove(output)

	print 'Lese Template ein'

	hTemplate = open(template_path, 'r')
	try:
		text = hTemplate.read()
	finally:
		hTemplate.close()

	print 'Do maple stuff'
	os.system(MAPLE_PATH + ' < ' + maple_path) 

	print 'Ersetzungen durchfuehren'
	i = 0
	for tmp_file in tmp_files:
		text_rep = u'${%d}$' % i
		htmp_file = open(tmp_file, 'r')

		try:
			text_input = htmp_file.read()
			text = text.replace(text_rep, text_input)
			
		finally:
			htmp_file.close()
		i = i + 1

	for rep in reps:
		text = text.replace(rep[0], rep[1])

	file_write_to = open(output, 'w')
	file_write_to.write(text)
	file_write_to.close()
	

	for filename in tmp_files:
		if os.path.exists(filename):
			os.remove(filename)


#	template_file = open(template_path, 'r')

#	try:
#		text = template_file.read()
#	finally:
#		template_file.close()

#	for tmp_file in tmp_files:











	
