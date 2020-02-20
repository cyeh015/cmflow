from geom_3dface_utils import *

import unittest
import os

# this file's location, use .join() to locate other files at the same directory
#     os.path.join(__location__, 'abc.xyz')
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

class TestFace(unittest.TestCase):
    """docstring for TestFace"""

    def test_load(self):
        fault = Face3D()
        fault.read(os.path.join(__location__, 'test_geom_3dface_utils_ex.ts'))
        self.assertEqual(len(fault.points),8286)
        self.assertEqual(len(fault.triangles),15857)

        self.assertEqual(fault.points['1'].name,'1')
        self.assertEqual(
            [p.name for p in fault.triangles[0].points],
            ['10', '9', '8'])
        self.assertEqual(
            [(p0.name, p1.name) for p0,p1 in fault.triangles[0].edges()],
            [('10', '9'), ('10', '8'), ('9', '8')])

class TestTriangles(unittest.TestCase):
    """docstring for TestTriangles"""

    def test_edges(self):
        t = Triangle([1,3,5])

        self.assertEqual(
            t.edges(),
            [(1,3), (1,5), (3,5)])

        self.assertEqual(
            t.edges([3]),
            [(1,5)])

        self.assertEqual(
            t.edges([(1,5)]),
            [(1,3), (3,5)])

class TestLineIntersectPlane(unittest.TestCase):
    """docstring for TestDataSeries"""

    def test_line_intersect_plane(self):
        # 1 pt, in seg
        self.assertEqual(
            line_intersect_plane((11.,5.,1.), (19.,25.,-1.), (10.,10.,0.), (0.,0.,1.),
                in_segment=False, check_on_plane=True),
            [15.0, 15.0, 0.0])
        # 1 pt, end point
        self.assertEqual(
            line_intersect_plane((11.,5.,1.), (19.,25.,0.), (10.,10.,0.), (0.,0.,1.),
                in_segment=True, check_on_plane=True),
            1)
        # 1 pt, the other end point
        self.assertEqual(
            line_intersect_plane((11.,5.,0.), (19.,25.,1.), (10.,10.,0.), (0.,0.,1.),
                in_segment=True, check_on_plane=True),
            0)
        # on plane, check
        self.assertEqual(
            line_intersect_plane((11.,5.,0.), (19.,25.,0.), (10.,10.,0.), (0.,0.,1.),
                in_segment=False, check_on_plane=True),
            (0,1))
        # on plane, no check
        self.assertEqual(
            line_intersect_plane((11.,5.,0.), (19.,25.,0.), (10.,10.,0.), (0.,0.,1.),
                in_segment=False, check_on_plane=False),
            None)
        # outside seg
        self.assertEqual(
            line_intersect_plane((11.,5.,1.), (19.,25.,10.), (10.,10.,0.), (0.,0.,1.),
                in_segment=True, check_on_plane=True),
            None)
        # outside seg, no extension
        self.assertEqual(
            line_intersect_plane((11.,5.,1.), (7.,-5.,2.), (10.,10.,0.), (0.,0.,1.),
                in_segment=False, check_on_plane=True),
            [15.0, 15.0, 0.0])
        # parallel
        self.assertEqual(
            line_intersect_plane((11.,5.,1.), (19.,25.,1.), (10.,10.,0.), (0.,0.,1.),
                in_segment=False, check_on_plane=True),
            None)

class TestFaultFullCase(unittest.TestCase):
    def test_heavy_real_case(self):
        """
        TODO: should setup lighter case and check solution, this barely checks
        if it runs ok.
        """
        fault = Face3D()
        fault.read(os.path.join(__location__, 'test_geom_3dface_utils_ex.ts'))
        print fault

        from mulgrids import mulgrid
        geo = mulgrid(os.path.join(__location__, 'test_geom_3dface_utils_geom.dat'))

        for lay in geo.layerlist:
            print '----- layer:', lay.name, lay.centre
            fault.set_cutting_plane((0.0,0.0,lay.centre), (0.0,0.0,1.0))
            line = [fault.search_line()]


if __name__ == '__main__':
    unittest.main(verbosity=2)


