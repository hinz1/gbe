#!/usr/bin/python
import sys
import argparse 
import os.path
import time

parser = argparse.ArgumentParser()
parser.add_argument('--gene', '-g', type=str, nargs='?', required=True, help="Use this Gene (like ATnG12345 or something).")
parser.add_argument('--basedir','-d', type=str, nargs='?', default='/mnt/work2/aschmidt/final_2018-12-13', help="Directory to work in, defaults to /mnt/work2/aschmidt/final_2018-12-13")
parser.add_argument('--out','-o', type=str, nargs='?', default='mutstat3_', help="Name of output file, defaults to mutstat3_ATnG12345.png")
parser.add_argument('--disp','-dp', action='store_true', default=False, help="Display the output instead of saving the image file? Defaults to false.")
parser.add_argument('--posperbar','-p', type=int, nargs='?', default=60, help="Aggregate how may positions per bar? Defaults to 60")
parser.add_argument('--width','-w', type=int, nargs='?', default=2400, help="Pixel width of figure. Defaults to 2400")
parser.add_argument('--aspect','-a', type=int, nargs='?', default=66, help="Height as percentage of width. Defaults to 66")
parser.add_argument('--cdsref','-r', type=str, nargs='?', default='Bostr', help="String to set CDS-reference, defaults to Bostr")
parser.add_argument('--debug','-dbg', action='store_true', default=False, help="Print debug output.")
parser.add_argument('--outfm','-of', type=str, nargs='?', default='svg', help="Name of output format, defaults to svg, could be e.g. png")
parser.add_argument('--cut5','-c5', type=int, nargs='?', default=0, help="Skip n bp upstream. Defaults to 0")
parser.add_argument('--cut3','-c3', type=int, nargs='?', default=0, help="Skip n bp downstream. Defaults to 0")
parser.add_argument('--summary','-s', action='store_true', default=False, help="Show summary plot.")
parser.add_argument('--sumsep','-sp', action='store_true', default=False, help="Show summary plot separately.")

args = parser.parse_args()

if args.sumsep == True:
	plts = 3
else:
	plts = 2
smplt = args.summary

g = args.gene
d = args.basedir
o = args.out
outfm = args.outfm
if o == 'mutstat3_':
	o += g + '.' + outfm
w = args.width
ht = float(args.aspect)/100
dpi = 100
dbg = args.debug
disp = args.disp
ppb = args.posperbar 
cut5 = args.cut5
cut3 = args.cut3

inpath = d + '/' + g + '/'
cdsnm = args.cdsref
aln = {}
cds = ''
if dbg:
	print(args)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

if os.path.isdir(inpath):
	if os.path.isfile(inpath + g + '.mafft.fa'):
		alnfile = inpath + g + '.mafft.fa'
	elif os.path.isfile('./' + g + '.mafft.fa'):
		alnfile = './' + g + '.mafft.fa'
	else:
		print('Could not find alignment: "' + g + '.mafft.fa' + '" or "' + inpath + g + '.mafft.fa' + '".')
		sys.exit()

	if os.path.isfile(alnfile):
		f = open(alnfile)
		a = f.read().splitlines()
		nseq = 0
		nseqa = 0
		nseqs = 0
		nseqf = 0
		fstlen = len(a[1])

		for sq in a:
			if sq[0] == '>':
				nseq += 1
				if '_a_' in sq:
					nseqa += 1
				elif '_s_' in sq:
					nseqs += 1
				elif '_fa_' in sq:
					nseqf += 1
				nm = sq[1:].strip()
				aln[nm] = ''
			else:
				aln[nm] += sq[cut5:fstlen-cut3]
				if cdsnm in nm:
					cds = sq[cut5:fstlen-cut3]
		maxlen = len(aln[nm])
		f.close()
	else:
		print('Could not read alignment: "' + alnfile + '".')
		sys.exit()
else:
	print('Could not find base dir: "' + inpath + '".')
	sys.exit()
pts = int(maxlen/ppb)
pos = 0
mutstata = [0] * (maxlen + ppb)
mutstats = [0] * (maxlen + ppb)
mutstatf = [0] * (maxlen + ppb)
mutstatges = [0] * (maxlen + ppb)
mutstatcn = [0] * (maxlen + ppb)
mutstatcountges = [0] * (maxlen + ppb)
mutstatcounta = [0] * (maxlen + ppb)
mutstatcounts = [0] * (maxlen + ppb)
mutstatcountf = [0] * (maxlen + ppb)
msctges = [0] * (maxlen + ppb + 2)
mscta = [0] * (maxlen + ppb + 2)
mscts = [0] * (maxlen + ppb + 2)
msctf = [0] * (maxlen + ppb + 2)

n = 0
for b in cds:
	b = b.lower()
	nmuta = 0
	nmuts = 0
	nmutf = 0
	cola = ''
	cols = ''
	colf = ''
	col = ''
	inex = 0.9
	if b == '-':
		inex = 0.0
	for s in aln:
		if cdsnm not in s:
			bb = aln[s][n].lower()
			#if b != bb:
			if '_a_' in s:
				if bb not in cola:
					nmuta += 1
					cola += bb
			elif '_s_' in s:
				if bb not in cols:
					nmuts += 1
					cols += bb
			elif '_fa_' in s:
				if bb not in colf:
					nmutf += 1
					colf += bb
			col += bb
	
	col = ''.join(set(col))
	#mutstata[n] = [inex, nmuta*100/nseqa, len(col)-1]
	#mutstats[n] = [inex, nmuts*100/nseqs, len(col)-1]
	#mutstatf[n] = [inex, nmutf*100/nseqf, len(col)-1]
	mutstata[n] = [inex, nmuta, len(col)-1]
	mutstats[n] = [inex, nmuts, len(col)-1]
	mutstatf[n] = [inex, nmutf, len(col)-1]
	mutstatges[n] = [inex, len(col), len(col)-1]

	mutstatcn[n] = mutstata[n][0]
	mutstatcountges[n] = mutstatges[n][1]
	mutstatcounta[n] = mutstata[n][1]
	mutstatcounts[n] = mutstats[n][1]
	mutstatcountf[n] = mutstatf[n][1]
	if dbg:
		if nmuta != nmuts or nmuta != nmutf or nmuts != nmutf:
			print(str(n) + ': ' + str(nmuta)+' '+str(nmuts)+' '+str(nmutf) +' ; '+ str(mutstatcounta[n]) + ' ' + str(mutstatcounts[n]) + ' ' + str(mutstatcountf[n]) + ' max: ' + str(maxlen))
	n += 1

n = 0
for msc in mutstatcounta:
	mscta[n] =   sum(mutstatcounta[n:n+ppb-1])/float(ppb)
	mscts[n+1] = sum(mutstatcounts[n:n+ppb-1])/float(ppb)+0.01
	msctf[n+2] = sum(mutstatcountf[n:n+ppb-1])/float(ppb)+0.02
	msctges[n] = sum(mutstatcountges[n:n+ppb-1])/float(ppb)
	n += 1
if dbg:
	print('a')
	print(mscta[100:500])
	print('s')
	print(mscts[100:500])
	print('fa')
	print(msctf[100:500])
	print('ges')
	print(msctges[100:500])

#fig = plt.figure(figsize=(w/dpi, w*0.66/dpi), dpi=dpi)
fig, axs = plt.subplots(plts, 1, figsize=(w/dpi, w*ht/dpi), dpi=dpi)
#axs[0].plot(mutstatcount)
axs[0].plot(msctf, color="green")
axs[0].plot(mscta, color="blue")
axs[0].plot(mscts, color="red")
if plts == 3:
	axs[1].plot(msctges, color="black")
	axs[2].plot(mutstatcn)
	axs[2].set_ylim(0, 1)
else: 
	axs[1].plot(mutstatcn)
	if smplt == True:
		axs[0].plot(msctges, color="black")
	axs[1].set_ylim(0, 1)
plt.suptitle('Mutations in ' + g) 
if disp == True:
	print("Showing plot.")
	plt.show()
else:
	print("Saving to " + o)
	plt.savefig(o, bbox_inches='tight', format=outfm)
