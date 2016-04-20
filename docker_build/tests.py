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

