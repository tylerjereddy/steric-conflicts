import unittest
import numpy as np
import scipy
import numpy.testing
import steric_assessment_general
import MDAnalysis

class TestVesicleStericViolations(unittest.TestCase):
    '''Verify properties of steric analysis using the DPPC vesicle example from the MARTINI website.'''

    @classmethod
    def setUpClass(cls):
        cls.u = MDAnalysis.Universe('dppc_vesicle.gro')
        cls.DPPC_selection = cls.u.select_atoms('resname DPPC')  
        cls.num_DPPC_residues = cls.DPPC_selection.n_residues #should be 877 DPPC residues
        cls.start_index = 1
        cls.end_index = cls.DPPC_selection.n_atoms
        cls.particles_per_residue = 12

    @classmethod
    def tearDownClass(cls):
        del cls.u
        del cls.DPPC_selection 
        del cls.num_DPPC_residues 
        del cls.start_index 
        del cls.end_index 
        del cls.particles_per_residue 

    def test_data_structure_steric_violations(self):
        '''Test for proper length of ordered_list_steric_violation_counts.'''
        list_DPPC_steric_violations = steric_assessment_general.main(start_index = self.start_index,end_index = self.end_index, coordinate_file = 'dppc_vesicle.gro',particles_per_residue = self.particles_per_residue, cutoff = 8.0)
        self.assertEqual(len(list_DPPC_steric_violations), self.num_DPPC_residues, "The number of DPPC steric violation counts should match the total number of DPPC residues in the vesicle system.")

    def test_cutoff_sensitivity(self):
        '''Test to ensure that using a larger cutoff leads to far more steric conflicts counted per residue.'''
        list_DPPC_steric_violations_small_cutoff = steric_assessment_general.main(start_index = self.start_index,end_index = self.end_index, coordinate_file = 'dppc_vesicle.gro',particles_per_residue = self.particles_per_residue, cutoff = 2.0)
        list_DPPC_steric_violations_large_cutoff = steric_assessment_general.main(start_index = self.start_index,end_index = self.end_index, coordinate_file = 'dppc_vesicle.gro',particles_per_residue = self.particles_per_residue, cutoff = 20.0)
        self.assertGreater(np.array(list_DPPC_steric_violations_large_cutoff).sum(), np.array(list_DPPC_steric_violations_small_cutoff).sum(), "The sum total per residue steric violations should be larger when using a larger cutoff.")

class TestIdenticalCoords(unittest.TestCase):
    '''Unit test(s) dealing with a simple coordinate file containing two identical copies of a DPPC molecule.'''

    def test_steric_violation_counts(self):
        '''Test the total number of steric violations for the identical residues (each should have 12 == number of particles in other residue.'''
        list_DPPC_steric_violations = steric_assessment_general.main(start_index = 1,end_index = 24, coordinate_file = 'dppc_simple_copies.gro',particles_per_residue = 12, cutoff = 2.0)
        self.assertEqual(len(list_DPPC_steric_violations), 2, "list_DPPC_steric_violations should have a length of 2 because there are only two residues to probe.")
        self.assertEqual(list_DPPC_steric_violations, [12,12], "Each DPPC molecule should have 12 steric violations, as both DPPC molecules in this test have the same coords.")

    def test_single_core_activity(self):
        '''Simple test to probe perform_portion_of_steric_analysis_on_individual_core() function.'''
        single_DPPC_list_steric_viols = steric_assessment_general.perform_portion_of_steric_analysis_on_individual_core('dppc_simple_copies.gro', 1, 12, 12, cutoff=2.0)
        self.assertEqual(single_DPPC_list_steric_viols, [12])

class TestAdjustArrays(unittest.TestCase):
    '''Test adjust_arrays_to_avoid_splaying() function, which has to be able to handle a variety of atom index data structures gracefully because different machines can have different numbers of cores.'''
    @classmethod
    def setUpClass(cls):
        #list index arrays based on two DPPC residues being probed across N cores
        cls.index_array = np.arange(1,25)
        cls.list_index_arrays_1_core = np.array_split(cls.index_array, 1)
        cls.list_index_arrays_2_cores = np.array_split(cls.index_array, 2)
        cls.list_index_arrays_3_cores = np.array_split(cls.index_array, 3)
        cls.list_index_arrays_7_cores = np.array_split(cls.index_array, 7)
        cls.list_index_arrays_12_cores = np.array_split(cls.index_array, 12)
        cls.list_index_arrays_19_cores = np.array_split(cls.index_array, 19)
        cls.list_index_arrays_36_cores = np.array_split(cls.index_array, 36)
        cls.list_index_arrays_70_cores = np.array_split(cls.index_array, 70)
 
        cls.list_index_arrays_expected_2_plus_cores = [np.arange(1,13), np.arange(13,25)] #there are only two residues, so never more than two arrays in the list irrespective of the number of cores
        cls.list_index_arrays_expected_1_core = [np.arange(1,25)] #for a single core, only a single array containing both residues can be distributed

    @classmethod
    def tearDownClass(cls):
        del cls.index_array 
        del cls.list_index_arrays_1_core 
        del cls.list_index_arrays_2_cores 
        del cls.list_index_arrays_3_cores 
        del cls.list_index_arrays_7_cores
        del cls.list_index_arrays_12_cores 
        del cls.list_index_arrays_19_cores
        del cls.list_index_arrays_36_cores 
        del cls.list_index_arrays_70_cores 
 
        del cls.list_index_arrays_expected_2_plus_cores 
        del cls.list_index_arrays_expected_1_core 

    def test_adjust_arrays_1_core(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_1_core, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_1_core):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 
    
    def test_adjust_arrays_2_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_2_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_3_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_3_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_7_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_7_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_12_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_12_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_19_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_19_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_36_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_36_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 

    def test_adjust_arrays_70_cores(self):
        list_index_arrays_obtained = steric_assessment_general.adjust_arrays_to_avoid_splaying(self.list_index_arrays_70_cores, 12)
        for actual_array, desired_array in zip(list_index_arrays_obtained, self.list_index_arrays_expected_2_plus_cores):
            np.testing.assert_allclose(actual_array, desired_array, rtol=1e-7) 
