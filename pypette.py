from __future__ import print_function
import sys, os, subprocess

###########
# LOGGING #
###########

def info(message):
	print(message, file=sys.stderr)


def error(message):
	print(message, file=sys.stderr)
	sys.exit(-1)


#########################
# FILESYSTEM MANAGEMENT #
#########################

def zopen(path, mode='r'):
	gzip_executable = 'gzip'
	#gzip_executable = 'pigz'
	
	if path[0] == '~':
		path = os.path.expanduser(path)
	
	if path == '-':
		if mode == 'r':
			return sys.stdin
		elif mode == 'w':
			return sys.stdout
	
	if path.lower().endswith('.gz'):
		if mode == 'r':
			return subprocess.Popen('gunzip -c %s' % path,
				stdout=subprocess.PIPE, shell=True).stdout
		elif mode == 'w':
			return subprocess.Popen('%s -c > %s' % (gzip_executable, path),
				stdin=subprocess.PIPE, shell=True).stdin
	else:
		return open(path, mode)


######################
# PROCESS MANAGEMENT #
######################

def shell(command, **kwargs):
	try:
		subprocess.check_call(command, shell=True, executable='/bin/bash',
			**kwargs)
	except subprocess.CalledProcessError as e:
		print('Process returned with error %d.' % e.returncode)
		return False
	return True

def shell_stdin(command):
	return subprocess.Popen(command, stdin=subprocess.PIPE, shell=True,
		executable='/bin/bash').stdin

def shell_stdout(command):
	return subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,
		executable='/bin/bash', bufsize=-1).stdout
