#!/bin/env python

import sys
import subprocess

proc = subprocess.Popen(['gnuplot', '-'], shell=True, stdin=subprocess.PIPE,)
proc.stdin.write('set term x11\n')
proc.stdin.write('set mouse\n')

doubletsFile = sys.argv[1]

with open(doubletsFile) as file:
    doublets = []
    for line in file:
        low,high,startTime,endTime,startScan,endScan = line.strip().split(',')
        doublets.append((low, high, float(startTime), float(endTime), startScan, endScan))

with open("promising.txt", "w") as outfile:
    index = 0
    maxDoublets = len(doublets)
    while index < maxDoublets:
        thisDoublet = doublets[index]
        low, high, startTime, endTime, startScan, endScan = thisDoublet
        msg1 = "doublet {}/{}: ".format(index,maxDoublets-1)
        msg2 = "%s %s spanning %.3f %.3f" % (low, high, startTime, endTime)
        print(msg1 + msg2)
        width = endTime - startTime
        lowX = startTime - (width * 3)
        highX = endTime + (width * 3)
        proc.stdin.write('set xrange [{}:{}]\n'.format(lowX, highX))
        proc.stdin.write('unset arrow\n')
        proc.stdin.write('set arrow from {0},graph 0 to {0},graph 1 nohead lw 0.5 lc rgb "grey"\n'.format(startTime))
        proc.stdin.write('set arrow from {0},graph 0 to {0},graph 1 nohead lw 0.5 lc rgb "grey"\n'.format(endTime))
        proc.stdin.write('plot "{0}" using 1:3 with lines ti "{0}", "{1}" using ($1):($3) with lines ti "{1}"\n'.format(low,high))
        input = raw_input("Enter for next, p for previous, q to quit, s to save: ").lower()
        if input == 'p':
            index -= 1
        elif input == 'q':
            proc.stdin.write('quit\n')
            break
        elif input == 's':
            outfile.write(','.join([low, high, str(startTime), str(endTime), str(startScan), str(endScan)]) + '\n')
            outfile.flush()
            index += 1
        else:
            index += 1
