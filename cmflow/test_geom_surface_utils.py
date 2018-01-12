import unittest
from geom_surface_utils import *

class TestSnap(unittest.TestCase):
    def setUp(self):
        self.sequence = [
            #  6      5      4      3     2       1        0
            100.00, 75.00, 50.00, 25.00, 0.00, -50.00, -100.00,
        ][::-1]
    def test_standard(self):
        self.assertEqual(25.0, snap(self.sequence, 27.0, None, None))
        self.assertEqual(25.0, snap(self.sequence, 27.0, 'down', None))
        self.assertEqual(50.0, snap(self.sequence, 27.0, 'up', None))
        self.assertEqual(50.0, snap(self.sequence, 27.0, 'down', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 27.0, 'down', 51.0))
        self.assertEqual(50.0, snap(self.sequence, 27.0, 'up', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 27.0, 'up', 51.0))
        self.assertEqual(50.0, snap(self.sequence, 38.0, None, None))
        self.assertEqual(25.0, snap(self.sequence, 38.0, 'down', None))
        self.assertEqual(50.0, snap(self.sequence, 38.0, 'up', None))
        self.assertEqual(50.0, snap(self.sequence, 38.0, 'down', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 38.0, 'down', 51.0))
        self.assertEqual(50.0, snap(self.sequence, 38.0, 'up', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 38.0, 'up', 51.0))
    def test_tricky(self):
        self.assertEqual(25.0, snap(self.sequence, 25.0, None, None))
        self.assertEqual(25.0, snap(self.sequence, 25.0, 'down', None))
        self.assertEqual(25.0, snap(self.sequence, 25.0, 'up', None))
        self.assertEqual(50.0, snap(self.sequence, 25.0, 'down', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 25.0, 'down', 51.0))
        self.assertEqual(50.0, snap(self.sequence, 25.0, 'up', 29.0))
        self.assertEqual(75.0, snap(self.sequence, 25.0, 'up', 51.0))
        self.assertEqual(100.0, snap(self.sequence, 100.0, None, None))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'down', None))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'up', None))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'down', 29.0))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'down', 51.0))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'up', 29.0))
        self.assertEqual(100.0, snap(self.sequence, 100.0, 'up', 51.0))
        self.assertEqual(-100.0, snap(self.sequence, -100.0, None, None))
        self.assertEqual(-100.0, snap(self.sequence, -100.0, 'down', None))
        self.assertEqual(-100.0, snap(self.sequence, -100.0, 'up', None))
        self.assertEqual(  50.0, snap(self.sequence, -100.0, 'down', 29.0))
        self.assertEqual(  75.0, snap(self.sequence, -100.0, 'down', 51.0))
        self.assertEqual(  50.0, snap(self.sequence, -100.0, 'up', 29.0))
        self.assertEqual(  75.0, snap(self.sequence, -100.0, 'up', 51.0))
    def test_cant_do(self):
        # if out of range
        self.assertRaises(Exception, snap, self.sequence, -105.0, None, None)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'down', None)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'up', None)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'down', 29.0)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'down', 51.0)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'up', 29.0)
        self.assertRaises(Exception, snap, self.sequence, -105.0, 'up', 51.0)
        self.assertRaises(Exception, snap, self.sequence, 105.0, None, None)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'down', None)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'up', None)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'down', 29.0)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'down', 51.0)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'up', 29.0)
        self.assertRaises(Exception, snap, self.sequence, 105.0, 'up', 51.0)

if __name__ == '__main__':
    unittest.main()
