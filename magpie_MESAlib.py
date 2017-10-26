import pandas as pd
import shutil
import os
import numpy as np
from subprocess import PIPE, Popen
_mesa_dir = os.environ.get('MESA_DIR')
_mesasdk_root = os.environ.get('MESASDK_ROOT')
_maxmultipars = 10000

class parameterGroup(object):
    """class that handles all parameters for a run 
    for browsing and editing inlists.
    
    Attributes:
        params (dict): used parameters in the run separated 
            by type ('controls', 'star_job', 'pgstar')
        defaults (dict): every parameter available in the 
            defaults files. Also split into types.
    """
    types = ['controls', 'star_job', 'pgstar']
    
    def __init__(self):
        self.params = {}
        self.defaults = {}
        for t in parameterGroup.types:
            self.params[t] = {}
            self.defaults[t] = {}

    def __getitem__(self, key):
        if key in parameterGroup.types:
            return self.params[key]
        else:
            return {}

    def readDefaults(self):
        """Read in all default parameters from the 
        star/defaults folder. star_job.defaults has the best 
        behaved syntax so docstrings should be correct.
        Controls and pgstar however, are as random as they come 
        so docstrings are seldom correct or available.
        """
        if _mesa_dir is None:
            print("MESA_DIR not set, cannot read defaults.")
            return
        globstr = _mesa_dir+"/star/defaults/"
        fnames = sorted([x for x in os.listdir(globstr) if ".defaults" in x])
        fpaths = [globstr+x for x in fnames]
        assert(len(fpaths)==3)
        self.setDefs(fpaths[0], type='controls')
        self.setDefs(fpaths[1], type='pgstar')
        self.setDefs(fpaths[2], type='star_job')
        for t in parameterGroup.types:
            self.mixPars(t)
            
    def setPars(self, inlist, type='controls'):
        """sets the parameters from an inlist in the 
        object. Reads in only the 'type' parameters so if 
        an inlist has multiple param types, run 'readInlist' 
        instead.
        
        Args:
            inlist (str): inlist to be read.
            type (str): type of parameter to read in 
                ('controls', 'star_job', 'pgstar').
        """
        cd = dict(getRawPars(inlist, type=type))
        self.params[type].update(cd)

    def setDefs(self, parfile, type='controls'):
        """reads in a .default file parameters and tries to 
        get docstrings for each parameter.
        
        Args:
            parfile (str): .default file to be read.
            type (str): type of parameter to read in 
                ('controls', 'star_job', 'pgstar').
        """
        n, v, d, err = getRawDefaults(parfile)
        cd = dict(zip(n, zip(v,d)))
        self.defaults[type].update(cd)

    def readInlist(self, inlist):
        """reads in all parameters from an inlist
        
        Args:
            inlist (str): inlist to be read.
        
        """
        self.setPars(inlist, type='controls')
        self.setPars(inlist, type='star_job')
        self.setPars(inlist, type='pgstar')
        
    def quickLook(self):
        """returns pandas dataframe with read inlist values for quick 
        editing of parameters.
        
        Returns:
            pandas.core.frame.DataFrame: tabulated params (all types).
        
        """
        frames = []
        for t in parameterGroup.types:
            df = pd.DataFrame()
            names, values = zip(*self.params[t].items())
            if not names:
                names, values = "None", 0
            df[t], df['value'] = names, values
            frames.append(df)
        return pd.concat(frames)
        
    def tabulate(self, type):
        """returns a pandas dataframe with parameters.
        if defaults have been read, returns every parameter available 
        with values from read inlist. This is a browsing mode.
        if there's no defaults, show currently used parameters.
        
        Args:
            type (str): type of parameters to show 
                ('controls', 'star_job', 'pgstar').
        """
        df = pd.DataFrame()
        docs = ""
        if not self.defaults[type]:
            names, values = zip(*self.params[type].items())
        else:
            self.mixPars(type)
            names, tuples = zip(*self.defaults[type].items())
            values, docs = zip(*tuples)
            docs = parseDocs(docs)
        if not names:
            names, values, docs = "None", 0, ""
        df[type], df['value'], df['docstr'] = names, values, docs
        return df.set_index(type)

    def mixPars(self, type='controls'):
        """sets the values from inlists into the default 
        dictionary.
        """
        if not self.defaults[type]:
            return
        else:
            for k, v in self[type].items():
                try:
                    dv, doc = self.defaults[type][k]
                    self.defaults[type][k] = (v, doc)
                except KeyError:
                    continue

    def readChanges(self, df):
        """Takes an edited pandas dataframe and reads 
        changed values.
        """
        for t in parameterGroup.types:
            try:
                cont = df[[t,'value']].dropna(0).set_index(t)
            except KeyError:
                if df.index.name==t:
                    cont = df
                else:
                    continue
            newdict = cont.T.to_dict('records')[0]
            for k, v in newdict.items():
                try:
                    dv, doc = self.defaults[t][k]
                    if v!=dv:
                        self.defaults[t][k] = (v, doc)
                        self[t][k] = v
                except KeyError:
                    continue

    def writeInlists(self, outpath, sjname="star_job",
                     ctname="controls", pgname="pgstar"):
        """Writes inlist files at 'outpath' using the object's params
        values.
        
        Args:
            outpath (str): path to work folder.
            sjname (str): star_job inlist suffix.
            ctname (str): controls inlist suffix.
            pgname (str): pgstar inlist suffix.
        """
        for i, t in enumerate(parameterGroup.types):
            keys = ["read_extra_{}_inlist1", "extra_{}_inlist1_name"]
            vals = ['.true.', "inlist_{}".format(t)]
            inlistd = dict(zip([k.format(t) for k in keys], vals))
            if not i:
                write2Inlist(inlistd, "&{}".format(t), outpath, "inlist", 
                             clobber=True)
            else:
                write2Inlist(inlistd, "&{}".format(t), outpath, "inlist")
            write2Inlist(self[t], "&{}".format(t), outpath, vals[1], 
                         clobber=True)


def parseDocs(dlist):
    """converts multiline docstrings into a single string."""
    nd = []
    for d in dlist:
        if isinstance(d, list):
            nd.append("".join(d))
        else:
            nd.append(d)
    return nd


def setupMESArun(destination, clobber=True):
    """copy the star/work folder from your MESA distro  to 
    a new folder to make a new run.
    
    Args:
        destination (str): folder path for new run.
        clobber (bool): overwrite the folder.
    
    """
    codesource = os.path.join(_mesa_dir, "star/work")
    if os.path.exists(destination):
        if clobber:
            shutil.rmtree(destination)
            shutil.copytree(codesource, destination)
    else:
        shutil.copytree(codesource, destination)


def getRawPars(inlist, type='controls'):
    """ gets a list of parameters from an inlist, 
    subject to type.
    Args:
        inlist (str): input inlist.
        type (str): header type for parameters.
    Returns:
        list: [(name, value), ...]
    """
    pars = []
    add = 0
    head = '&{}'.format(type)
    with open(inlist, 'r') as parf:
        for line in parf:
            if len(line)<=1:
                continue
            else:
                if head in line:
                    add = 1
                    continue
                elif not add:
                    continue
                else:
                    if line[0]=='/':
                        break
                    if '=' in line and '!' not in line:
                        l = line.strip().split('=')
                        pars.append((l[0].strip(),l[-1].strip()))
    if not pars:
        return []
    else:
        return pars


def getRawDefaults(defaultfile):
    """ returns list of names, defaults and descriptions found for each 
    parameter in a .defaults file. if syntax is incorrect the parameter 
    won't have a description but the name will still be returned for use.
    Fills 'err' with the lines around the wrong syntax parameters.
    Args:
        defaultfile (str): .default file.
    
    Returns:
       list: names, list: values, list: docstrings, list: err.
    """
    pars, desc, multi = [], [], []
    name, val = "", ""
    doc, setvars = False, False
    pnames, err = [], []
    flush = False
    with open(defaultfile, 'r') as par:
        for line in par:
            if len(line.strip())>1:
                if line.strip()[0]!= '!':
                    prel =line.split('=')[0].strip('9012345678 .):({}#!\n')
                    pnames.append(prel.lower())
            if len(line.strip())<1:
                continue
            if flush and not '!###' in line.split()[0]:
                err.append(line.strip('\n'))
                continue
            if '!###' in line.split()[0]:
                flush = False
                if setvars:
                    #this means there's not enough values for the defined 
                    # variables, so reset everything
                    name, val = "", ""
                    desc = []
                    multi = []
                    setvars = False
                if not name:
                    # split fixes comentaries right beside the DEFINITION 
                    # of the variables. lower fixes Msun vs msun differences.
                    # no consistency whatsoever :C
                    prel = line.strip('9012345678  .):({}#!\n').lower()
                    name = prel.split()[0]
                    multi.append(name)
                    namenum = plainGenerator(_maxmultipars)
                else:
                    prel = line.strip('9012345678  .):({}#!\n').lower()
                    multi.append(prel.split()[0])
                doc = True
                continue
            if doc:
                if line.strip()[0] == '!':
                    desc.append(line.strip(' !'))
                else:
                    setvars = True
                    # controls.defaults has double '=' on some of its params...
                    avoidequals = line.strip().split('=')
                    cname, val = avoidequals[:2]
                    cname = cname.strip('9012345678  .):(').lower()
                    if len(multi)<=1:
                        try:
                            m = multi.index(cname)
                        except ValueError:
                            flush = True
                            err.append(line.strip('\n'))
                            continue
                        pars.append((multi[m], val, desc))
                        name, val = "", ""
                        desc = []
                        multi = []
                        doc = False
                        setvars = False
                    else:
                        n = namenum.next()
                        try:
                            m = multi.index(cname)
                        except ValueError:
                            flush = True
                            err.append(line.strip('\n'))
                            continue
                        if n==len(multi)-1:
                            pars.append((multi[m], val, u"".join(desc)))
                            name, val = "", ""
                            desc = []
                            multi = []
                            doc = False
                            setvars = False
                        else:
                            pars.append((multi[m], val, "".join(desc)))
            else:
                continue
    tags = [a[0] for a in pars]
    for st in pnames:
        if st not in tags:
            pars.append((st, "None", "None"))
    n, v, d = zip(*pars)
    return n, v, d, err


def plainGenerator(length):
    """simple integer counting generator.
    
    Args:
        length (int): max count for generator.
        
    Yields:
        int: counter from 0 to 'length'.
    
    """
    i = 0
    while i<length:
        yield i
        i += 1


def fortParse(arg):
    """returns a parsed variable from a parameter (bool,
    str, or number)
    
    Args:
        arg (str): parameter value
        
    Returns:
        str: decorated argument for fortran parsing.
    
    """
    try:
        val = np.float(arg.replace('d','E'))
        return arg
    except ValueError:
        if arg=='.true.':
            return arg
        elif arg=='.false.':
            return arg
        else:
            return '"{}"'.format(arg.strip('"\' '))


def write2Inlist(parameters, header, outpath, fname, clobber=False):
    """Write paramater ditionary to file, appending for
    clobber=False (default)
    
    Args:
        parameters (dict): dictionary of parameters.
        header (str): type of parameters.
        outpath (str): MESA work folder.
        fname (str): inlist filename.
        clobber (bool): rewrite file or append to file.
    
    """
    if clobber:
        opt = 'w'
    else:
        opt = 'a'
    with open("/".join([outpath, fname]), opt) as o:
        o.write("\n{}\n\n".format(header))
        for k, v in parameters.items():
            o.write("      {} = {}\n".format(k, fortParse(v)))
        o.write("\n/\n")


def compileMESA(outpath, startsdk=True):
    """calls ./mk at outpath, compiling the run
    
    Args:
        outpath (str): MESA work folder path.
        startsdk (bool): source the mesasdk before compiling.
    
    """
    if _mesa_dir is None:
        print("MESA_DIR not set. Cannot compile. Returning.")
        return 1
    if startsdk:
        if _mesasdk_root is None:
            print("MESASDK_ROOT not set. Skipping sdk init.")
            comm = 'cd {} && ./mk'.format(outpath)
        else:
            init = '{}/bin/mesasdk_init.sh'.format(_mesasdk_root)
            comm = 'source {} && cd {} && ./mk'.format(init, outpath)
    print('Excecuting: {}'.format(comm))
    p = Popen(['/bin/bash'], stdin=PIPE, stdout=PIPE)
    out, err = p.communicate(input=comm.encode())
    print(out.decode())
    exitcode = p.returncode
    return exitcode