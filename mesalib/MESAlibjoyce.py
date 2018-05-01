#!/usr/bin/env python
#try:
import numpy as np
import codecs, re
import subprocess, os
#import h5py
#import pygadgetreader as pyg
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
#import healpy as hp
import random as rand

################################################################
#
# the purpose of this module is MESA HANDLING
#
################################################################

#except:
#    print 'Missing module!\nThe following are required: '
#    print 'healpy\n'
#   exit(0)


def plotter(xmaj, xmin, ymaj, ymin,**kwargs):
    xf = kwargs.get('xf', '%1.1f')
    yf = kwargs.get('yf', '%1.1f')

    majorLocator_x  = MultipleLocator(xmaj)     # I want a major tick every "number"
    majorFormatter_x = FormatStrFormatter(xf)#('%1.1f')     # label these (the major ones) with a 1.2f format string
    minorLocator_x  = MultipleLocator(xmin)     # I want a minor tick every "number"

    majorLocator_y  = MultipleLocator(ymaj)     # now for the y axis...
    majorFormatter_y = FormatStrFormatter(yf)#('%1.1f')     # 
    minorLocator_y  = MultipleLocator(ymin)     #

    fig, ax = plt.subplots()

    ax.xaxis.set_major_locator(majorLocator_x)
    ax.xaxis.set_major_formatter(majorFormatter_x)
    ax.xaxis.set_minor_locator(minorLocator_x)

    ax.yaxis.set_major_locator(majorLocator_y)
    ax.yaxis.set_major_formatter(majorFormatter_y)
    ax.yaxis.set_minor_locator(minorLocator_y)

    return fig,ax

###
# ideally need dictionary or some hash table containing map of i for p[i] to physical quantity name 
### even better than this: read in all the column titles from the MESA output in one shot, shove them in a list,
# iterate over that list using enumerate or whatever, then BUILD a dictionary out of names AS LOADED FROM THE THING and assign an index to that 
# simply by iteration


def submit_MESA_job(mesa_work_path,qsub_script):
    mesa_work_path=str(mesa_work_path)
    os.chdir(mesa_work_path)
    #	subprocess.call('cd '+mesa_work_path, shell=True)
    print "entered"
    subprocess.call('pwd', shell=True)
    compile_MESA()
    run_MESA()
    return 

def compile_MESA():
	try:
		subprocess.call('./mk', shell=True)
	except:
		print 'MESA star compilation failed'
	return

def execute_MESA():
	subprocess.call('cd '+mesa_work_path, shell=True)
	try:
		subprocess.call('./rn', shell=True)
	except:
		print 'MESA star execution failed'	
	return


def update_MESA_inlist_value(MESA_inlist, field, value):
    inlist_dict=grab_fields(MESA_inlist)
    old_value=inlist_dict.get(field)
    new_value=str(value)
    oldstr=r'{}\s*=\s*{}'.format(field, old_value)
    newstr='{} = {}'.format(field, value)
    #print("oldstr: {}\nnewstr: {}".format(oldstr, newstr))
    f = codecs.open(MESA_inlist)
    contents = f.read()
    newcontents=re.sub(oldstr, newstr, contents) 
    f.close()
    outf=open(MESA_inlist,"w")
    print >>outf, newcontents
    outf.close()
    return 

def grab_fields(MESA_inlist):
    inf=open(MESA_inlist,'r')
    inlist_dict={}
    for line in inf:
        if line and "=" in line:
            try:
                p=line.split("=")
                try: 
                    p[1]=p[1].split('\n')[0]
                except:
                    pass
                if "!" in p[0]:
                    del p    
                inlist_dict[str(p[0]).strip()]=str(p[1]).strip()
            except:
                pass 
    return inlist_dict


#########################################################################
#
# MESA output functions
#
#########################################################################
def strip_MESA_header(in_filename, out_filename, *args, **kwargs):

    ### add hook to prevent double removes

	#[0]-- returns the file object
	#[1]-- returns the string/name of the reformatted text file
	n = int(kwargs.get('n', 5))
	num_delete=int(n) #the number of line to be read and deleted
	outfn=out_filename

	with open(in_filename) as f:
	    mylist = f.read().splitlines()
	newlist = mylist[:]
	#os.remove("bigFile.txt")

	thefile = open(outfn, 'w')
	del mylist[:num_delete]
	for item in mylist:
		thefile.write("%s\n" % item)
	return  thefile, outfn

def get_MESA_output_fields(filename):
    inf=open(filename,'r')
    for line in inf:
        #print "for line in inf"
        # not robust but whatever for now
        #try:
            #print "in try catch"
        #print line
        if "zone" in line and "num_zones" not in line:
            #try:
            #print line
            p=line.split()
            #except:
            #    strip_MESA_header(filename,filename)
            #    p=line.split()
        #except:
        else:
            #print "history"
            if 'luminosity' in line:                    
                #try:
                p=line.split()
                #except:
                #    strip_MESA_header(filename,filename)
                #    p=line.split()
    phys_dict={}
    try:		 	
        for i,v in enumerate(p):
            phys_dict[str(v)]=i#-1 	#wow fucking damn if this was the problem 
            #yes, it was.
    except UnboundLocalError:
        print "problem with "+str(filename)+" format"
        sys.exit(0)  
    return phys_dict


def get_columns(filename,keyname_list):
    phys_dict=get_MESA_output_fields(filename)
    column_dict={}
    for i in range(len(keyname_list)):
        try:
            keyname=str(keyname_list[i])
            column_dict[ keyname ]=(get_key(filename, keyname)[0])
        except:
            print show_allowed_MESA_keywords(readfile)
    return column_dict


def get_key(filename,keyname):
    inf=open(filename,'r')
    phys_dict=get_MESA_output_fields(filename)
    indx=phys_dict.get(str(keyname))
    param=[]
    for line in inf:
      #avoid 'zone' instead of '#'
      if line and 'zone' not in line:
               #print line
          try: 
              p=line.split()
              param.append(p[indx])
          except:
              pass
    return param, type(param)

def get_quantity(readfile,keyname):
    keyname=str(keyname)
    keyname_list=get_MESA_output_fields(readfile).keys()
    column_dict=get_columns(readfile,keyname_list)
    quantity=np.array(column_dict.get(keyname))[3:]# NO WRONG .astype(float)
    #magic 3 to eliminate column numbers being interpreted as data in the profile file
    return np.array(quantity).astype(float)


def show_allowed_MESA_keywords(readfile):
    fstr=""
    for i in get_MESA_output_fields(readfile).keys():
        fstr=fstr+str(i) +'\n'
    print fstr
    return fstr#get_MESA_output_fields(readfile).keys()




def transfer_inlist(template_inlist, target_inlist):
    with open(template_inlist, 'r') as content_file:
        contents = content_file.read()

    outf=open(target_inlist,"w")
    print >> outf, contents
    outf.close()

    return



def generate_basic_inlist(mass, age, metallicity, m2g_path, mesa_dir, inlist_path,output_model_name):
    ## I don't know how paths work
    outf=open(inlist_path, "w")

    print >> outf, "&star_job"
    print >> outf, ""
    print >> outf, "  mesa_dir = '"+ str(mesa_dir) +"'"
    print >> outf, "  history_columns_file='"+m2g_path+"/data/history_columns_testsuite.list'"
    print >> outf, "  profile_columns_file='"+m2g_path+"/data/profile_columns_testsuite.list'"
    print >> outf, ""
    print >> outf, "  load_saved_model = .false."
    print >> outf, "  !create_pre_main_sequence_model = .false. !!need this to be able to specify y?"
    print >> outf, ""
    print >> outf, "  save_model_when_terminate = .false."
    print >> outf, ""
    print >> outf, "  write_profile_when_terminate = .true."
    print >> outf, "  filename_for_profile_when_terminate = 'profile_"+str(output_model_name)+".data'"
    print >> outf, ""
    print >> outf, "  pgstar_flag = .false."
    print >> outf, "  save_pgstar_files_when_terminate = .false."
    print >> outf, ""
    print >> outf, " / !end of star_job namelist"
    print >> outf, ""
    print >> outf, "&controls"
    print >> outf, ""
    print >> outf, "  star_history_name = 'history_"+str(output_model_name)+".data'"
    print >> outf, "" 
    print >> outf, "  !max_number_backups = 0"
    print >> outf, "  max_number_retries = 5"
    print >> outf, "  !max_model_number = 500"
    print >> outf, ""
    print >> outf, "  initial_mass ="+str(mass)+"d0"
    print >> outf, "  !initial_y = 0.27d0 !solar"
    print >> outf, "  initial_z = "+str(metallicity)+"d0 !solar"
    print >> outf, ""
    print >> outf, "  terminal_show_age_in_years = .true."
    print >> outf, "  max_age = "+str(age)+"d9 ! age of Sun in Gyr"
    print >> outf, ""
    print >> outf, ""
    print >> outf, "  !photo_interval = 50"
    print >> outf, "  !profile_interval = 50"
    print >> outf, "  history_interval = 1"
    print >> outf, "  !terminal_interval = 100"
    print >> outf, "  !write_header_frequency = 100"
    print >> outf, ""
    print >> outf, "/ ! end of controls namelist"
    print >> outf, ""
    outf.close()
    return 
