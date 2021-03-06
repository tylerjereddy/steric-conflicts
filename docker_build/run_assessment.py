'''Author: Tyler Reddy
   Purpose: To perform a multicore steric assessment using the multiprocessing-based general-purpose module.'''

import MDAnalysis
import numpy as np
import cPickle as pickle
import steric_assessment_general 
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

parser = argparse.ArgumentParser()
parser.add_argument("-start_index", type=int)
parser.add_argument("-end_index", type=int)
parser.add_argument("-coord_filepath", type=str)
parser.add_argument("-particles_per_residue", type=int)
parser.add_argument("-cutoff", type=float, help="cutoff (A)")
parser.add_argument("-pickle_filename", type=str, default='steric_viols.p', nargs='?')
parser.add_argument("-plot_filename", type=str, default='steric_histogram.png', nargs='?')
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
matplotlib.pyplot.xticks(rotation=0)
unique_contact_counts, num_residues_with_those_counts = np.unique(array_steric_violations, return_counts=True)
percent_denominator = array_steric_violations.size / 100.0
ax.bar(unique_contact_counts, num_residues_with_those_counts / percent_denominator,facecolor = 'green',alpha=0.75,width=0.7, align='center') #percent histogram
ax.set_xlim(-1, unique_contact_counts.max() + 1)
ax.set_xticks(np.arange(-1,unique_contact_counts.max() + 2))
ax.set_xlabel('# of contacts within {cutoff} $\AA$'.format(cutoff=args.cutoff))
ax.set_ylabel('% of residues')
fig.set_size_inches(6,6)
fig.subplots_adjust(bottom=0.2, left=0.2)
fig.savefig(args.plot_filename,dpi=300)
