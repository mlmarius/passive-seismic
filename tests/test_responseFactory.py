import pytest

import os, sys
from os.path import join

from InstrumentResponse import ResponseFactory
import numpy as np

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
ASDF = join(PASSIVE, 'ASDFdatabase')
TESTDATA = join(ASDF, 'testdata')

def test_ResponseFromPAZ():
    rf = ResponseFactory()

    # Create a Response Object from poles and zeros and associated variables that define
    # a system response.
    poles = [-59.4313+0j, -22.7121+27.1065j, -22.7121-27.1065j, -0.0048004+0j, -0.073199+0j]
    zeros = [1+0.33j,2+0.45j]
    rf.CreateFromPAZ('trillium', 'LAPLACE (RADIANS/SECOND)',
                 86083,
                 0.02,
                 1935,
                 0.02,
                 poles,
                 zeros)

    r = rf.getResponse('trillium')
    assert(np.allclose(np.array(r.get_paz().poles), np.array(poles), 1e-5))
    assert(np.allclose(np.array(r.get_paz().zeros), np.array(zeros), 1e-5))
#end func

def test_ResponseFromStationXML():
    rf = ResponseFactory()

    # Create a Response Object a StationXML file, derived from a RESP file
    rf.CreateFromStationXML('le3dlite', join(TESTDATA, 'LE3dLite.xml'))

    r = rf.getResponse('le3dlite')
    assert(np.allclose(np.array(r.get_paz().poles), np.array([-4.444+4.444j, -4.444-4.444j, -1.083+0.j]), 1e-5))
    assert(np.allclose(np.array(r.get_paz().zeros), np.array([ 0.+0.j,  0.+0.j,  0.+0.j]), 1e-5))
#end func