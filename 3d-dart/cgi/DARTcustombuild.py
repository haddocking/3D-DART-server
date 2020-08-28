#!/usr/bin/env python2.7

USAGE = """
==========================================================================================

Author:               Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
                      for Biomolecular Research, Utrecht university, The Netherlands
Copyright (C):        2012 (DART project)
DART version:         1.3  (17-01-2012)
DART system module:   dart.py
Module function:      DART html form wrapper script
Plugin dependencies:  Standard python and CGI module (included in Python)

==========================================================================================
"""

"""import modules"""
import cgi, sys, os
import cgitb; cgitb.enable(display=1, logdir='/var/www/html/3DDART/error')

"""Setting variables"""
DARTDIR = (os.getcwd()).replace('cgi','server')

if DARTDIR not in sys.path:
	sys.path.append(DARTDIR)

sys.path.append(DARTDIR+"/system/")

from system.DARTserver import WebServer

#Get the full http server path to service
if os.environ.has_key('SCRIPT_URI'):
	HTTPDIR = os.path.dirname(os.environ['SCRIPT_URI'].replace('/cgi',''))
else:
	HTTPDIR = '/3DDART/html/'


class DartCustomBuild(object):	
	def __init__(self):
		# Get values from webform
		self.webform = cgi.FieldStorage(keep_blank_values=True)
		# Initialize default python dictionary
		self.python_dict = {}
		
		self.__cgiToPythonDict()
		self.__processWebform()
		
	def __checkDictKeys(self, keys, default=''):
		"""Checks the dictionary if it has the keys. If not add them with value None"""
		for n in keys:
			if not self.python_dict.has_key(n):
				self.python_dict[n] = default

	def __cgiToPythonDict(self):
		"""Retrieve all keys and value pairs from the CGI FieldStorage class. This makes
		   postprocessing easier as the FieldStorage class does not allow keys to be added.
		"""
		for keys in self.webform:
			try:
				# For single value entries for one parameter
				self.python_dict[keys] = (self.webform[keys].value).strip()
			except:
				# For multiple value entries for one parameters
				self.python_dict[keys] = (",".join(self.webform.getvalue(keys, "")))
	
	def __checkFilePath(self,path):
		"""Check to see if file path contains backslashes. This can happen if user uploads file from Internet Explore.
		   What a suprise, MS always does it in a different way. If backslahes are found, convert them to forewards 
		   slashes and feed the string to os.path.basename and return result
 		"""
		tmp = path.replace('\\','/')
		if tmp:
			return os.path.basename(tmp)
		else:
			# Multiple values are formated to comma seperated list
			return path
			
	def __processWebform(self):
		"""Make some custom build specific checks for the variables past to the workflow.xml file"""
		# If a parameter or pdb file is uploaded, build from this one and blank the sequence entry field. 
		# If a sequence is provided check its base-pair count. More than 200 base-pairs are not 
		# allowed to prevent server crash.
		seq = False
		self.__checkDictKeys(['1.upload'])
		if len(self.python_dict['2.sequence']):
		    if (len(self.python_dict['2.sequence']) * int(self.python_dict['2.repeat'])) > 1000:
		        self.__returnError("Too many base-pairs (%s): A nucleic acid sequence of more than 1000 base-pairs is not supported" 
		                            % (len(self.python_dict['2.sequence']) * int(self.python_dict['2.repeat'])))
		    seq = True		    
		if len(self.python_dict['1.uploadpar']):
			self.python_dict['2.sequence'] = ''
			self.python_dict['2.repeat'] = '1'
			self.python_dict['1.upload'] = {'file':self.python_dict['1.uploadpar'],
											'name':self.__checkFilePath(self.webform['1.uploadpar'].filename)}
			seq = True
		elif len(self.python_dict['1.uploadpdb']): 
			self.python_dict['2.sequence'] = ''
			self.python_dict['1.upload'] = {'file':self.python_dict['1.uploadpdb'],
											'name':self.__checkFilePath(self.webform['1.uploadpdb'].filename)}
			# Disable first BuildNucleicAcids plugin
			self.python_dict['2.useplugin'] = 'False'
			# Have the third X3DNAanalyze plugin use the uploaded PDB
			self.python_dict['3.inputfrom'] = '1.0'
			
			#Check if the uploaded file is a .zip file. If so, change workflow to ensemblebuild
			if os.path.splitext(self.python_dict['1.upload']['name'])[1] == '.zip':
				self.python_dict['NAensemblebuild.xml'] = 'submit'
				del self.python_dict['NAcustombuild.xml']
			seq = True	
		
		# If no sequence and no file upload report error
		if not seq:
			self.__returnError("Empty fields: Please provide either a nucleic acid sequence a\
								custom base-pair(step) parameter file or a PDB structure file to model DNA from")
		
		# NA type is set once and applies to NAbuild job2 and job5
		self.__checkDictKeys(['2.type','5.type'])
		self.python_dict['5.type'] = self.python_dict['2.type']
	
		# Once Calladin-Drew plot in block 1 or block 2 representation always that representation
		self.__checkDictKeys(['2.block1','5.block1','2.block2','5.block2'],default='False')
		self.python_dict['5.block1'] = self.python_dict['2.block1']
		self.python_dict['5.block2'] = self.python_dict['2.block2']
		
		# Reference base-pair for generation is same as for analysis		
		self.__checkDictKeys(['4.refbp','7.refbp'])
		self.python_dict['7.refbp'] = self.python_dict['4.refbp']
	
		# If a restraints file need to be made than pdb2haddock has to be true
		self.__checkDictKeys(['9.useplugin','8.pdb2haddock'])
		if self.python_dict['9.useplugin'] == 'True':
			self.python_dict['8.pdb2haddock'] = 'True'
		else:
			self.python_dict['9.useplugin'] = 'False'
	
		# Test to see if the user wants to model anything other than canonical DNA.
		# If so then the fiber structure is used for the PDBeditor stage and not the 3DNA rebuild structure
		Model_test = ['4.anglerange','4.orientrange','4.bpstep','4.bp']
		Wform_haskey = True
		NothingToModel = True
	
		for test in Model_test:
			if not self.python_dict.has_key(test):
				Wform_haskey = False
	
		if Wform_haskey == True:
			for test in Model_test:
				if not self.python_dict[test] == '' or not len(self.python_dict[test]) == 5:
					NothingToModel = False
			if NothingToModel == True:
				self.__checkDictKeys(['8.inputfrom'])
				self.python_dict['8.inputfrom'] = '2.0'
		
		#If something needs to be modelled, extract start and end basepair from angelzone and set angle range
		# and orient range
		if not NothingToModel:
			#Set min/max bend angle
			if self.python_dict.has_key('4.anglerange'):
				if self.python_dict['4.global'] == 'True':
					if '-' in self.python_dict['4.anglerange']:
						splitted = self.python_dict['4.anglerange'].split('-')
						if len(splitted) == 2:
							self.python_dict['4.minangle'] = splitted[0].strip()
							self.python_dict['4.maxangle'] = splitted[1].strip()
					else:
						self.python_dict['4.minangle'] = self.python_dict['4.anglerange']
				else:
					self.python_dict['4.minangle'] = self.python_dict['4.anglerange'].replace('-',',')
			del self.python_dict['4.anglerange']
			
			#Set min/max bend orientation
			if self.python_dict.has_key('4.orientrange'):
				if self.python_dict['4.global'] == 'True':
					if '-' in self.python_dict['4.orientrange']:
						splitted = self.python_dict['4.orientrange'].split('-')
						if len(splitted) == 2:
							self.python_dict['4.minorient'] = splitted[0].strip()
							self.python_dict['4.maxorient'] = splitted[1].strip()
					else:
						self.python_dict['4.minorient'] = self.python_dict['4.orientrange']
				else:
					self.python_dict['4.minorient'] = self.python_dict['4.orientrange'].replace('-',',')
			del self.python_dict['4.orientrange']
			
			#Set start and end bp for bending
			if self.python_dict.has_key('4.anglezone'):
				if '-' in self.python_dict['4.anglezone']:
					splitted = self.python_dict['4.anglezone'].split('-')
					if len(splitted) == 2:
						self.python_dict['4.startbp'] = splitted[0].strip()
						self.python_dict['4.endbp'] = splitted[1].strip()
					else:
						self.python_dict['4.startbp'] = splitted[0].strip()
				else:
					self.python_dict['4.startbp'] = self.python_dict['4.anglezone']
			del self.python_dict['4.anglezone']
		
	def __returnError(self,error):
		"""If an traceable error occured present it to the user and then exit"""
		print('Location: %s/DARTerror.html?err=%s\n' % (HTTPDIR, error))
		sys.exit()
		
	def createEnv(self):
		"""Gather information about the user for log statistics"""
		environment = []
		for x in ['HTTP_USER_AGENT','REMOTE_ADDR','REMOTE_HOST','REMOTE_PORT','HTTP_REFERER']:
			if os.environ.has_key(x):
				environment.append("%s: %s\n" % (x, os.environ[x]))
		return environment			

	def dartDownload(self, filepath=None):
		"""Present the download path of the data to the user"""
		print('Status: 303 See other')
		print('Location: %s/DARTresult.html?loc=%s\n' % (HTTPDIR, filepath))


if __name__ == '__main__':

	custom = DartCustomBuild()
	
	# Send form to webserver
	server = WebServer(DARTDIR, remote_env=custom.createEnv())
	
	# Run DART from in webserver mode
	filepath = server.RunDART(custom.python_dict)
	
	# Collect download path and present to user
	custom.dartDownload(filepath)
