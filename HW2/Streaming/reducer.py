#!/usr/bin/env python

from operator import itemgetter
import sys
import csv

current_line = None
current_count = 0
line = None

for line in sys.stdin:
    line = line.strip()
    count = int(1)

    if current_line == line:
        current_count += count
    else:
        if current_line:
            x, y  = line.split('_')
            x_low = str(float(x)-0.1)
            y_low = str(float(y)-0.1)
            print '%s\t%s\t%s\t%s\t%s' % (x_low, x,  y_low, y, str(current_count))
        current_count = count
        current_line = line

if current_line == line:
    x, y  = line.split('_')
    x_low = str(float(x)-0.1)
    y_low = str(float(y)-0.1)
    print '%s\t%s\t%s\t%s\t%s' % (x_low, x,  y_low, y, str(current_count))









