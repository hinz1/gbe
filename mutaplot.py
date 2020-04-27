#!/usr/bin/python
import sys
import argparse 
import os.path
import math

parser = argparse.ArgumentParser()
parser.add_argument('--gene', '-g', type=str, nargs='?', required=True, help="Use this Gene (like ATnG12345 or something).")
#parser.add_argument('--max', '-mx', type=int, nargs='?', default=5000, help="Maximum length. defaults to 5000")
parser.add_argument('--maf', '-m', type=str, nargs='?', default="allcov50.maf", help="Use this maf file, defaults to /mnt/work2/aschmidt/formafffinal/allcov50.maf")
parser.add_argument('--width', '-w', type=int, nargs='?', default=2400, help="Width of image, defaults to 2400.")
parser.add_argument('--basedir','-d', type=str, nargs='?', default='/mnt/work2/aschmidt/formafffinal', help="Directory to work in, defaults to /mnt/work2/aschmidt/formafffinal")
parser.add_argument('--out','-o', type=str, nargs='?', default='mutaplot', help="Name of output file, defaults to mutaplot_ATnG12345.png")
parser.add_argument('--disp','-dp', action='store_true', default=False, help="Display the output instead of saving the image file? Defaults to false.")
#parser.add_argument('--nseq','-n', type=int, nargs='?', default=145, help="")
parser.add_argument('--posperbar','-p', type=int, nargs='?', default=60, help="Aggregate how may positions per bar? Defaults to 60")
parser.add_argument('--height','-ht', nargs='?', default=0.66, help="Fraction of image height of width, defaults to 0.66")
parser.add_argument('--outfm','-of', type=str, nargs='?', default='svg', help="Name of output format, defaults to svg, could be e.g. png")
parser.add_argument('--filter','-f', type=str, nargs='?', default='', help="String to use as search filter, e.g. '_a_', defaults to '_'")

args = parser.parse_args()

d = args.basedir
w = args.width
out = args.out
outfm = args.outfm
mafin = args.maf
gene = args.gene
flt = args.filter
#maxlen = args.max
#nseq = args.nseq
disp = args.disp
agg = args.posperbar
ht = float(args.height)
dpi=100
#print(args.height)
#print(ht)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

if os.path.isdir(d):
	if os.path.isfile(d + '/' + gene + '.mafft.cds.fa'):
		alnpath = d + '/' + gene + '.mafft.cds.fa'
	elif os.path.isfile(gene + '.mafft.cds.fa'):
		alnpath = gene + '.mafft.cds.fa'
	else:
		print('Could not find alignment: "' + gene + '.mafft.cds.fa' + '" or "' + d + '/' + gene + '.mafft.cds.fa' + '".')
		sys.exit()

	if os.path.isfile(d + '/' + gene + '.mafft.cds.fa'):
		f = open(alnpath)
		aln = f.read().splitlines()
		fstln = aln[1]
		maxlen = len(fstln)
		nseq = 0
		for sq in aln:
			if sq[0] == '>' and flt in sq:
				nseq += 1
		f.close()
	else:
		print('Could not read alignment: "' + d + '/' + gene + '.mafft.cds.fa' + '".')
		sys.exit()
else:
	print('Could not find base dir: "' + d + '".')
	sys.exit()

bars = int(maxlen/agg) #wieviele Balken
bw = 0.9
#print(fstln)
#print('Gene: ' + str(gene))
#print('mxlen: ' + str(maxlen))
#print('MAF: ' + str(mafin))
#print('n seq: ' + str(nseq))
#print('Bars: ' + str(bars) + ' a ' + str(bw))
x = range(1, bars+1)
ynst = [0] * bars
ynse = [0] * bars
yfs = [0] * bars
yif = [0] * bars
ysil = [0] * bars
ynstp = [0] * bars
ynsep = [0] * bars
yfsp = [0] * bars
yifp = [0] * bars
ysilp = [0] * bars

if os.path.isfile(d + '/' + mafin):
	mafinpath = d + '/' + mafin
elif os.path.isfile(mafin):
	mafinpath = mafin
else:
	print('Could not find MAF file: "' + mafin + '" or "' + d + '/' + mafin + '".')
	sys.exit()

f = open(mafinpath, 'r')
maf = f.read().splitlines()
f.close()

for mut in maf:
	#print(mut)
	mut = mut.split("\t")
	if len(mut) > 3:
		if mut[1] == gene and flt in mut[0]: 
			if mut[2] == 'Nonstop_Mutation': 
				ynst[int(round(int(mut[3])*bars/maxlen))] += 1
				ynstp[int(round(int(mut[3])*bars/maxlen))] = math.ceil(ynst[int(round(int(mut[3])*bars/maxlen))] * 100 / float(agg*nseq))
				# ynstp: agg*nseq=max. moegl. Anz. Mut., wieviel Prozent gefunden?
			elif mut[2] == 'Nonsense_Mutation': 
				ynse[int(round(int(mut[3])*bars/maxlen))] += 1
				ynsep[int(round(int(mut[3])*bars/maxlen))] = math.ceil(ynse[int(round(int(mut[3])*bars/maxlen))] * 100 / float(agg*nseq))
			elif mut[2] == 'Frame_Shift_Del' or mut[2] == 'Frame_Shift_Ins':
				frm = int(round(int(mut[3])*bars/maxlen))
				ntl = int(round(int(mut[4])*bars/maxlen))
				for i in range(frm, ntl+1):
					yfs[i] += 1
					yfsp[i] = math.ceil(yfs[i] * 100 / float(agg*nseq))
			elif mut[2] == 'In_Frame_Del' or mut[2] == 'In_Frame_Ins':
				frm = int(round(int(mut[3])*bars/maxlen))
				ntl = int(round(int(mut[4])*bars/maxlen))
				for i in range(frm, ntl+1):
					yif[i] += 1
					yifp[i] = math.ceil(yif[i] * 100 / float(agg*nseq))
			elif mut[2] == 'Silent': 
				ysil[int(round(int(mut[3])*bars/maxlen))] += 1
				ysilp[int(round(int(mut[3])*bars/maxlen))] = math.ceil(ysil[int(round(int(mut[3])*bars/maxlen))] * 100 / float(agg*nseq))
		
#plt.plot(x, ynst, 'rs')
#plt.plot(x, ynse, 'bs')

#print(x)
#print(ysil)
#print('Nonstop')
#print(ynstp)
#print('In frame')
#print(yifp)
#print('Nonsense')
#print(ynsep)
#print(ynse)
#print('Frameshift')
#print(yfsp)
#print(yfs)
#print(bars)
#print(disp)

fig = plt.figure(figsize=(w/dpi, w*ht/dpi), dpi=dpi)

bt = [0] * bars
lsil = plt.bar(x,ysilp,width=bw,bottom=bt,color='blue', edgecolor='black')
bt = [a+b for a,b in zip(bt,ysilp)]
lnst = plt.bar(x,ynstp,width=bw,bottom=bt,color='lime', edgecolor='black')
bt = [a+b for a,b in zip(bt,ynstp)]
lif  = plt.bar(x,yifp,width=bw,bottom=bt,color='darkorange', edgecolor='black')
bt = [a+b for a,b in zip(bt,yifp)]
lnse = plt.bar(x,ynsep,width=bw,bottom=bt,color='orangered', edgecolor='black')
bt = [a+b for a,b in zip(bt,ynsep)]
lfs  = plt.bar(x,yfsp,width=bw,bottom=bt,color='blueviolet', edgecolor='black')


plt.title('Gene ' + gene + ': percent mutations per site compared to reference, ' + str(agg) + ' sites per bar')
xt = [i*agg for i in x]
j = 0
k = 0
for i in xt:
	if j < 50:
		xt[k] = ''
	else:
		xt[k] = str(i)
		j = 0
	j+=1
	k+=1

plt.xticks(x, xt)
plt.legend((lsil[0], lnst[0], lif[0], lnse[0], lfs[0]), ('Silent', 'Non-stop', 'In-frame InDel', 'Non-sense', 'Frameshift'))
if disp == True:
	print("Showing plot.")
	plt.show()
else:
	outpath = d + '/' + out + '_' + gene + flt + '.' + outfm
	if os.access(d, os.W_OK):
		print("Saving to " + outpath)
		plt.savefig(outpath, bbox_inches='tight', format=outfm)
	else:
		print('Could not write output file "' + outpath + '"')
