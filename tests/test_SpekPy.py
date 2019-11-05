# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################

import spekpy as sp
from spekpy.IO import write_json_to_disk, read_json_from_disk
import pytest
import numpy as np
import sys, os

# Directory for test files and data
testdir='tests'

# Function to convert a class to a dictionary
def dict_from_class(cls):
    class A(object):
        pass
    _excluded_keys = set(A.__dict__.keys())
    dict_out = dict((key, value) for (key, value) in cls.__dict__.items() if key not in _excluded_keys)
    return dict_out

# Function to compare two dictionaries, element by element
def compare_dicts(dict1, dict2, atol=1e-8):
    errors=[]
    for item in dict1:
        ans1 = dict1[item]
        ans2 = dict2[item]
        if type(ans1) == np.float64:
            if not np.isclose(ans1, ans2, atol = atol):
                errors.append(item+': ref='+str(ans1)+' calc='+str(ans2))
        elif type(ans1) == np.ndarray:
            rmsd=np.sqrt(np.mean((ans1 - ans2)**2))
            if not np.isclose(rmsd, 0.0, atol=atol):
                errors.append(item+': rmsd(ref-calc)='+str(rmsd))
    return errors

# Function to generates reference values for the get_std_results() method of the Spek class
def generate_values_for_get_std_results(physics_choice):
    s=spek_ref(physics_choice)
    std_results_obj = s.get_std_results()
    dict1 = dict_from_class(std_results_obj)
    # Create a dictionary of results in a form that can be handled by json
    dict2 = {}
    for item in dict1:
        ans = dict1[item]
        if type(ans) == np.float64:
            dict2[item] = ans
        elif type(ans) == np.ndarray:
            dict2[item] = ans.tolist()
    # Save reference resullts using json
    write_json_to_disk(dict2,physics_choice + '_values_for_get_std_results.json')
    return

# Function to generates reference values for all methods of interest for the Spek Class
def generate_values(physics_choice):
    generate_values_for_get_std_results(physics_choice)
    return

# Reference spectrum definition
def spek_ref(physics_choice="default"):
    s = sp.Spek(th=8,physics=physics_choice).filter('Al',1.0)
    return s

# Code to generated new reference values if both (a) execute this script directly and (b) "generate" is first argument 
if __name__=="__main__":
    if sys.argv[1]=="generate":
        generate_values("default") # Generate values for physics = "default" mode
        generate_values("legacy")  # Generate values for physics = "legacy" mode    

# Test functions #
##################

# Test function for export_spectum() and load_from_file() methods of Spek class
def test_export_spectrum_and_load_from_file():
    s = spek_ref("default")
    spek_name = 'Test_export_spectrum.spk'
    # Reference values
    std_results_obj_1=s.get_std_results()
    s.export_spectrum(spek_name, comment='A test spectrum export')
    delim=';'
    e = sp.Spek.load_from_file(spek_name, delim)
    # Values after export and import
    std_results_obj_2=e.get_std_results()
    dict1=dict_from_class(std_results_obj_1)
    dict2=dict_from_class(std_results_obj_2)
    # Create a list of errors
    errors=compare_dicts(dict1,dict2,atol=1e-4)
    assert errors==[]


# Test function for get_std_results() method of Spek class in "legacy" mode
def test_get_std_results_legacy():
    physics_choice = "legacy"
    s = spek_ref(physics_choice)
    # Read in reference values
    dict1=read_json_from_disk(os.path.join(testdir,physics_choice+'_values_for_get_std_results.json'))
    std_results_obj=s.get_std_results()
    dict2=dict_from_class(std_results_obj)
    errors=compare_dicts(dict1,dict2)
    contents=dir(std_results_obj)   
    assert errors==[]

# Test function for get_std_results() method of Spek class in "default" mode
def test_get_std_results_default():
    physics_choice = "default"
    s = spek_ref(physics_choice)
    # Read in reference values
    dict1=read_json_from_disk(os.path.join(testdir,physics_choice+'_values_for_get_std_results.json'))
    std_results_obj=s.get_std_results()
    dict2=dict_from_class(std_results_obj)
    errors=compare_dicts(dict1,dict2)
    contents=dir(std_results_obj) 
    assert errors==[]

# Test function for save_state() and load_state() methods of the Spek class
def test_save_load_remove_state():
    s=spek_ref()
    std_results_obj_1=s.get_std_results()
    state_name='My spectrum state'
    s.save_state(state_name)
    t=sp.Spek.load_state(state_name)
    std_results_obj_2=t.get_std_results()
    sp.Spek.remove_state(state_name)
    dict1=dict_from_class(std_results_obj_1)
    dict2=dict_from_class(std_results_obj_2)
    errors=compare_dicts(dict1,dict2)
    assert errors==[]

# Test function to compare 1st HVLs of NIST reference spectra states (the "NIST_NIST" set) to values in the SpekPy paper
def test_nist_states_against_publication():
    paper_hvls=read_json_from_disk(os.path.join(testdir,'default_values_for_nist_hvls_in_paper.json'))
    errors = []
    for item in paper_hvls:
        hvlpaper=paper_hvls[item]
        s=sp.Spek.load_state(item)
        hvl=s.get_hvl1()
        if not np.isclose(hvlpaper,hvl,atol=1e-2):
            errors.append(item+': ref='+str(hvlpaper)+' calc='+str(hvl))
    assert errors == []

# Test function to check that the show_matls() method of the Spek class works with all keyword options
# Only tests that no exceptions are generated
def test_show_matls():
    sp.Spek.show_matls(matl_name='Water, Liquid')
    sp.Spek.show_matls(matl_group='ICRP')
    sp.Spek.show_matls(matl_group='ICRU')
    sp.Spek.show_matls(matl_dir='usr')
    sp.Spek.show_matls(matl_dir='def')
    sp.Spek.show_matls()

# Test function to check that the show_states() method of the Spek class works with all keyword options
# Only tests that no exceptions are generated
def test_show_states():
    sp.Spek.show_states(state_dir='usr')
    sp.Spek.show_states(state_dir='def')
    sp.Spek.show_states()

# Test function for make_matl() and remove_matl() methods of the Spek class
def test_make_matl():
    s=spek_ref()
    material_name = 'Wood, White Oak (test)'
    comment = 'PIET-43741-TM-963 PNNL-15870 Rev. 1 (http://www.pnnl.gov/main/publications/external/technical_reports/pnnl-15870rev1.pdf)'
    material_density = 0.77
    material_composition = [(1, 0.059642), (6, 0.497018), (7, 0.004970), (8, 0.427435), (12, 0.001988), (16, 0.004970), (19, 0.001988), (20, 0.001988)]
    sp.Spek.make_matl(matl_name=material_name, matl_density=material_density, wt_matl_comp=material_composition, matl_comment=comment)
    s.get_hvl1(matl=material_name)
    another_material = 'Water (body temperature) (test)'
    sp.Spek.make_matl(matl_name=another_material, matl_density=0.992, chemical_formula='H2O')
    s.get_hvl1(matl=another_material)
    sp.Spek.remove_matl(matl_name=material_name)
    sp.Spek.remove_matl(matl_name=another_material)

# Test function for whether teh set() method of the Spek class behaves as expected
def test_set():
    ref_spek='NIST_NIST_L80'
    s=spek_ref()
    s.set(kvp=80,th=21)
    s.filter('Air',1000).filter('Be',1.0).filter('Al',0.45)
    hvl=s.get_hvl1()
    k=s.get_kerma()
    t=sp.Spek.load_state(ref_spek)
    hvlref=t.get_hvl1()
    assert np.isclose(hvlref,hvl), 'set() produces inconsistent results'
    ref_kerma=100
    s.set(ref_kerma=ref_kerma,x=10,y=20)
    ktest=s.get_kerma()
    assert np.isclose(ktest,ref_kerma), 'Reference kerma normalization failed'
    ref_flu=100
    s.set(ref_flu=ref_flu,x=10,y=20)
    ftest=s.get_flu()
    assert np.isclose(ftest,ref_flu), 'Reference fluence normalization failed'
    s.set(x=0,y=0)
    kafter=s.get_kerma()
    assert np.isclose(k,kafter), 'Normalization not undone by subsequent call to set()'

# Test function to test filter() and multi_filter() methods of the Spek class are consistent
# and that removal of filters (negative values) can undo the addition of a filter
def test_filtering():
    s=spek_ref()
    t=spek_ref()
    std_results_obj_0=s.get_std_results()
    s.filter('Al',2.5).filter('Cu',0.25)
    std_results_obj_1 = s.get_std_results()
    t.multi_filter([('Al',2.5),('Cu',0.25)])
    std_results_obj_2 = t.get_std_results()
    dict1=dict_from_class(std_results_obj_1)
    dict2=dict_from_class(std_results_obj_2)
    errors=compare_dicts(dict1,dict2)
    assert errors==[], 'filter() gives different results to multi_filter()'
    t.multi_filter([('Al',-2.5),('Cu',-0.25)])
    std_results_obj_3 = t.get_std_results()
    dict0=dict_from_class(std_results_obj_0)
    dict3=dict_from_class(std_results_obj_3)
    errors=compare_dicts(dict0,dict3)
    assert errors==[], 'Removing filtration fails'

# Test function to test that the get_spectrum() method of the Spek class gives correct outputs
def test_get_spectrum():
    s=spek_ref()
    k0=s.get_k()
    spk0=s.get_spk()
    k, spk = s.get_spectrum()
    assert np.all(k==k0), 'get_spectrum() does not return identical energy bins to get_k()'
    assert np.all(spk==spk0), 'get_spectrum() does not return identical fluences to get_spk()'
    ke, spke = s.get_spectrum(edges='edges')
    assert np.size(k0)*2 == np.size(ke), 'get_spectrum(edges="edges") does not produce correct output'

# Test function to test that the get_matl() methods gives consistent and accurate results
def test_get_matl():
    s=spek_ref()
    hvl1=s.get_hvl1(matl='Cu')
    hvl2=s.get_hvl2(matl='Cu')
    t1=s.get_matl(matl='Cu',frac=0.5)
    assert np.isclose(hvl1,t1), 'get_matl() failed in "frac" mode'
    t2=s.get_matl(matl='Cu',hvl_matl='Cu',hvl=hvl2)
    assert np.isclose(t1,t2), '"frac" and "hvl" modes of get_matl() not consistent'

# Test function to test the clone method of the Spek class
def test_clone():
    s=spek_ref()
    hvl1=s.get_hvl1()
    t=sp.Spek.clone(s)
    t.filter('Al',hvl1)
    hvl2=t.get_hvl1()
    assert np.isclose(s.get_hvl2(), t.get_hvl1()), 'clone() does not behave as expected'
