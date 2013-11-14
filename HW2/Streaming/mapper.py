#!/usr/bin/env python

import sys

# write the function to round the value
from math import ceil, floor
def num_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)


# input comes from STDIN (standard input)

for line in sys.stdin:
    line = line.strip()
    # split the line into (group, value)
    number_float = map(float, line.split())
    r = []
    for i in range(0,len(number_float)):
        # write the results to STDOUT (standard output);
        # what we output here will be the input for the
        # Reduce step, i.e. the input for reducer.py
        r.append(num_round(number_float[i],1,ceil)) # round the value to the first decimal
    x = str(r[0])
    y = str(r[1])
    print '%s' % (x+'_'+y)

