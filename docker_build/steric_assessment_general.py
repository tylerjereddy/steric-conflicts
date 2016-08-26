'''
Author: Tyler Reddy
Purpose of module: Demonstrate utility of Python multiprocessing module for analysis of a single coordinate file with MDAnalysis. This module is also a useful general piece of code for analyzing a system that may have 'questionable contacts--' steric issues with atoms potentially too close together for any number of reasons (perhaps the most common use case would be after the construction of a complex system).
'''

import MDAnalysis
import multiprocessing
import time 
import sys
import numpy as np
import scipy
import cPickle as pickle

#We normally parse coordinate files with MDAnalysis as Universe objects, but these contain open files and cannot easily be passed between cores, so instead each core will launch the same analysis function with a different set of arguments, load a copy of the universe object, perform the analysis on the subset of the coordinates specified in the arguments, and finally return the results to the parent process for aggregation with results from other cores

def perform_portion_of_steric_analysis_on_individual_core(coord_file, species_start_index, species_end_index, particles_per_residue,cutoff=2.0):
    '''Each core will individually load the universe object (as it cannot be passed between them) and produce distance array data for a subset of the index range of interest.
    The arguments:
    coord_file : is a string containing the path to the single coordinate file being parsed
    species_start_index : is the atom index for the first particle to be parsed by THIS core
    species_end_index : is the atom index for the last particle to be parsed by THIS core
    particles_per_residue : is the number of particles in a given residue (i.e., 13 for coarse-grained POPC in MARTINI) so that the code can stride through and assess contacts around a single residue at a time    cutoff : the farthest a particle can be from a given residue to be counted as a proximal contact (measured in A)
    *Note: this code assumes that the indices of a given species are contiguous in the coordinate file.'''
    #use the multiprocessing module to print the name of the current process, so that we can monitor progress at the terminal:
    print multiprocessing.current_process().name, 'Starting' 
    #as usual with MDAnalysis, we produce a Universe object in order to do most of our work, and this will be done on each core individually as there's no easy way to pass Universe objects (which contain open files) between cores:
    input_universe = MDAnalysis.Universe(coord_file) #again, this will happen on every core that runs this function
    #we'll want to count the number of atoms that are within some cutoff of a given residue, and store all the values calculated on this core in the following list:
    dict_steric_violation_counts = {}
    #now we stride through the range of indices for this particular residue (i.e, POPC or DPPC, etc.) and count the number of atoms / particles that fall within the specified cutoff:
    current_species_index = species_start_index
    while current_species_index < species_end_index:
        current_species_residue_start_index = current_species_index #first atom in current residue
        current_species_residue_end_index = current_species_index + (particles_per_residue - 1) #last atom in current residue
        atom_selection_within_cutoff_angstroms = input_universe.select_atoms('around {cutoff} (bynum {start_index}:{end_index})'.format(start_index = current_species_residue_start_index, end_index = current_species_residue_end_index, cutoff = cutoff)) #simply select all atoms within cutoff A of the current residue
        number_of_atoms_within_cutoff = atom_selection_within_cutoff_angstroms.n_atoms #since my plan is to make a histogram, I'm really just interested in the number of violations for a given residue
        if number_of_atoms_within_cutoff > 0:
            print 'Violating residue start index:', current_species_residue_start_index ,'end index:', current_species_residue_end_index
        print multiprocessing.current_process().name, 'number of atoms within cutoff',number_of_atoms_within_cutoff #to monitor on a per-process basis at the terminal
        dict_steric_violation_counts[current_species_index] = number_of_atoms_within_cutoff
        current_species_index += particles_per_residue #stride forward to the starting index of the next residue
    print multiprocessing.current_process().name, 'Finishing' #to monitor end of this process on terminal
    return dict_steric_violation_counts

def adjust_arrays_to_avoid_splaying(start_index,end_index,particles_per_residue, available_cores = None):
    '''Adjust the arrays of indices sent to each core to avoid the splaying of residues across cores. After other code roughly balances out the number of indices to be sent to each core, this function takes the list of index arrays destined for each core and performs an adjustment such that the returned list contains index arrays which only have whole residues--NO splaying of residues between arrays/cores can be tolerated for sensible results. Naturally, it would be much simpler to deal with residue numbers rather than indices in simple cases, but index numbers are often more reliable in extremely large systems or in other coordinate files where there may be non-unique residue numbers.'''
    #I need to split the topology indices into contiguous chunks of residues (NO partial residues) so that different index ranges can be parsed by different cores:
    if available_cores == None: #detect automatically
        available_cores = multiprocessing.cpu_count() #determine number of CPUs available
    else: #specify available cores -- mostly for testing purposes
        available_cores = available_cores
    #start with a single array containing the full range of indices:
    index_array = np.arange(start_index,end_index + 1) #see numpy documentation
    #produce a list of arrays to be passed to the cores, with each array as evenly balanced as possible in terms of the number of indices it contains (use numpy array_split):
    num_indices = index_array.size
    num_residues = int(num_indices / float(particles_per_residue))

    if num_residues <= available_cores:
        available_cores = num_residues

    list_index_arrays = np.array_split(index_array,available_cores)
    #**note that an even (or uneven) workload balance between cores does NOT guarantee that residues are not splayed across cores; for this, I'll use the function defined above--which will shuffle the indices assigned to each core until each core is assigned an index range that contains only WHOLE residues (more important priority than a precise load balance):
    current_element = 0
    for index_array in list_index_arrays: #each index_array is destined for a core (later in the code)
        while index_array.size % particles_per_residue != 0: #don't exit the while loop until you have an integer number of residues accounted for in the index array
            index_array = np.concatenate((index_array,np.reshape(index_array[-1] + 1,(1,)))) #add the next index value
            #each time I add a value to the end of the above array I'll have to remove the first element of the subsequent array to avoid duplicating an index in my overall analysis:
            list_index_arrays[current_element + 1] = list_index_arrays[current_element + 1][1:]
        #it might seem like the above could fail on the last element (last array) in the list, but by definition that array should always be divisible by the particles_per_residue when we get to it because the overall number of indices is also divisible by the number of particles in a residue, so the remaining indices (last element / array) should not require modification
        list_index_arrays[current_element] = index_array #assign the new array, which only contains whole residues, to the list of arrays
        current_element += 1 #increment the array I'm working on after the current array is divisible by the particles per residue
    return list_index_arrays

#I'm now going to write a main control function for the module, so that when it is imported for a specific-use case on any arbitrary number of cores, the code will attempt to gracefully distribute the workload over the various cores:
def main(start_index,end_index,coordinate_file, particles_per_residue, cutoff):
    '''Main control function of the process. The arguments are as described for the per-core function above, except that the start and end indices are the overall start and end indices whereas the previous function receives an index range which is a subset of the index range specified here (this housekeeping work is dealt with by the code in this function).'''
    list_index_arrays_to_distribute_to_cores = adjust_arrays_to_avoid_splaying(start_index = start_index, end_index = end_index, particles_per_residue = particles_per_residue)
    #**want to be absolutely certain that each core receives a number of indices that is divisible by particles_per_residue:
    for index_array_for_core in list_index_arrays_to_distribute_to_cores:
        assert index_array_for_core.size % particles_per_residue == 0, "index arrays are splaying residues across cores"

    #finally, can start doing some of the multiprocessing heavy-lifting:
    pool = multiprocessing.Pool() #activate a pool of worker processes equal to the number of cores available on the system
    cumulative_parent_dict_steric_violation_counts = {} #this is where ALL per-residue steric violation counts will end up (in the parent process, progressively incorporating results from child processes)
    #note that the order of the results from children -- > parent is not known ahead of time in this module, but we could tag the results from each child in a dictionary, etc., if we had a workflow that required this

    def log_result_pool(dict_from_this_process): 
        '''This function will be called after each child process exits so that the overall dict of steric violations in the parent is grown by the dict of values produced in a given child.'''
        cumulative_parent_dict_steric_violation_counts.update(dict_from_this_process)

    #now, iterate through the list of index arrays and hand the steric assessment tasks off to the various cores:
    for index_array in list_index_arrays_to_distribute_to_cores:
        starting_index = index_array[0]
        ending_index = index_array[-1]
        #the apply_async function calls the function we wrote earlier for use by individual cores, and specifies which indices to use on a given core; the callback argument allows us to use the above logging function to dump the data back to the parent for aggregation purposes
        pool.apply_async(perform_portion_of_steric_analysis_on_individual_core, args = (coordinate_file, starting_index,ending_index,particles_per_residue,cutoff),callback = log_result_pool)
    #the next two methods basically ensure that the parent process waits for all the child processes to complete (see the multiprocessing docs for details)
    pool.close()
    pool.join()
    # generate a list of steric violations, sorted by the topological index of the residues in the coordinate file
    cumulative_parent_list_steric_violation_counts = []
    for index in sorted(cumulative_parent_dict_steric_violation_counts):
        cumulative_parent_list_steric_violation_counts.append(cumulative_parent_dict_steric_violation_counts[index])
    return cumulative_parent_list_steric_violation_counts #so if you use this code within another module, you'll just get the overall list of steric violations per residue
