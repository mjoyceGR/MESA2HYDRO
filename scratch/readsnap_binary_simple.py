#!/usr/bin/env python
import numpy as np
#import readsnap_binary_meridith as rbm
import io_lib as rbm


# snapshot filename
f=open('../data/snapshot_000',"r")

# particle type
ptype=1

header=rbm.load_gadget_binary_header(f)
attribute_dictionary=rbm.load_gadget_binary_particledat(f, header, ptype, skip_bh=0)

## shows which keyword to use to extract a specific physical property 
print attribute_dictionary.keys()

## example:
v=attribute_dictionary['Velocities']
positions=attribute_dictionary['Coordinates']

