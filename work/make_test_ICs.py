#!/usr/bin/env python
import numpy as np
import subprocess

subprocess.call("./run.py mainsequence.cfg", shell=True)
subprocess.call("./run.py redgiant.cfg", shell=True)
subprocess.call("./run.py OB.cfg", shell=True)
subprocess.call("./run.py AGB.cfg", shell=True)
subprocess.call("./run.py wd3.cfg", shell=True)

