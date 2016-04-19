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
