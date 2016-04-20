'''Author: Tyler Reddy
   Purpose: To perform a multicore steric assessment using the multiprocessing-based general-purpose module.'''

import MDAnalysis
import numpy as np
import cPickle as pickle
import steric_assessment_general 
import argparse
import matplotlib
import matplotlib.pyplot
matplotlib.use('Agg')

parser = argparse.ArgumentParser()
parser.add_argument("start_index", type=int)
parser.add_argument("end_index", type=int)
parser.add_argument("coord_filepath", type=str)
parser.add_argument("particles_per_residue", type=int)
parser.add_argument("cutoff", type=float, help="cutoff (A)")
parser.add_argument("pickle_filename", type=str, default='steric_viols.p')
parser.add_argument("plot_filename", type=str, default='steric_histogram.png')
args = parser.parse_args()

#produce an overall list of per-residue steric violations using the general multicore code in another module and convert to numpy array:
array_steric_violations = np.array(steric_assessment_general.main(start_index = args.start_index,
                                                                       end_index = args.end_index,
                                                                       coordinate_file = args.coord_filepath,
                                                                       particles_per_residue = args.particles_per_residue,
                                                                       cutoff = args.cutoff)) 

#because this may take quite a while, pickle the array to be safe:
pickle.dump(array_steric_violations,open(args.pickle_filename,'wb'))

#now plot the data (default plot -- the user can use the pickled data to customize their own):
fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(111)
matplotlib.pyplot.xticks(rotation=90)
histogram, bins = np.histogram(array_steric_violations,bins=20)
bincenters = 0.5 * (bins[1:] + bins[:-1]) #adjust to center them 
percent_denominator = array_steric_violations.size / 100.0
ax.bar(bincenters, histogram / percent_denominator,facecolor = 'green',alpha=0.75,width=2) #percent histogram
ax.set_xlim(-1,95)
ax.set_xlabel('# of contacts within {cutoff} $\AA$'.format(cutoff=args.cutoff))
ax.set_ylabel('% of residues')
fig.set_size_inches(4,6)
fig.savefig(args.plot_filename,dpi=300)
