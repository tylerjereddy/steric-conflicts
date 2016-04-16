'''Author: Tyler Reddy
   Purpose: To perform a multicore steric assessment using a standard test file along with the multiprocessing-based general-purpose module built above.'''

#import some useful modules:
import MDAnalysis, time, numpy, scipy, sys
import cPickle as pickle
import steric_assessment_general #you may need to ensure that the module is in your path with sys.path.append('path/to/module') for example

start_time = time.time() #start time in seconds for benchmarking purposes

#produce an overall list of steric violations using the general multicore code in another module and convert to numpy array:
array_DPPC_steric_violations = numpy.array(steric_assessment_general.main(start_index = 1,end_index = 10524,coordinate_file = 'dppc_vesicle.gro',particles_per_residue = 12, cutoff = 8.0)) #note that I've exaggerated the cutoff size a bit as the standard file doesn't likely have any major steric violations at, say, 2.0 A

#because this may take quite a while, pickle the array to be safe (see pickle docs for details):
pickle.dump(array_DPPC_steric_violations,open('DPPC_steric_viols.p','wb'))

#when simply adjusting the plot, can load directly from pickle with the above code commented (no need to re-run analysis code unless you want to change parameters):
#array_DPPC_steric_violations = pickle.load(open('DPPC_steric_viols.p','rb'))

#now plot the data (see matplotlib docs):
import matplotlib, matplotlib.pyplot

fig=matplotlib.pyplot.figure()
ax = fig.add_subplot(111)
matplotlib.pyplot.xticks(rotation=90)
histogram,bins = numpy.histogram(array_DPPC_steric_violations,bins=20)
bincenters = 0.5*(bins[1:]+bins[:-1]) #adjust to center them 
percent_denominator = array_DPPC_steric_violations.size / 100.0
ax.bar(bincenters,histogram / percent_denominator,facecolor = 'green',alpha=0.75,width=2) #percent histogram
ax.set_xlim(-1,95)
ax.set_xlabel('# of contacts within 8.0 $\AA$')
ax.set_ylabel('% of DPPC residues')
fig.set_size_inches(4,6)
fig.savefig('DPPC_steric_histogram.png',dpi=300)

print 'Steric assessment code completed in',time.time() - start_time, 'seconds'
