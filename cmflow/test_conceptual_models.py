from conceptual_models import *

import os
import unittest

class TestBMStats(unittest.TestCase):
    def test_new_save_load(self):
        geo1 = mulgrid().rectangular(
            [10.] * 10,
            [10.] * 10,
            [10.] * 10,
            convention=0,
            atmos_type=2, # no atm blocks
            origin=[0,0,0],
            justify='r',
            case=None,
            chars=ascii_lowercase)
        self.assertEqual(len(geo1.block_name_list), 1000)

        # new/empty bms requires geo
        with self.assertRaises(BMStatsError):
            bms = BMStats()

        # new/empty bms
        bms = BMStats(geo=geo1)
        self.assertEqual(bms.geo, geo1)
        self.assertEqual(bms.stats.shape, (1000, 0)) # np.array (num_blocks, num_zones)
        self.assertEqual(bms.zones, [])

        # save bms
        geo1.write('_bms1.dat')
        bms.save('_bms1.json')
        self.assertTrue(os.path.isfile('_bms1.json'))
        self.assertTrue(os.path.isfile('_bms1.dat'))
        self.assertTrue(os.path.isfile('_bms1.npy'))

        # load bms, with pre-loaded geo
        bms2 = BMStats('_bms1.json', geo=geo1)
        self.assertEqual(bms2.geo, geo1)
        self.assertEqual(bms2.stats.shape, (1000, 0))

        # load bms
        bms3 = BMStats('_bms1.json')
        self.assertEqual(bms3.geo.block_name_list, geo1.block_name_list)
        self.assertEqual(bms3.stats.shape, (1000, 0))
        self.assertEqual(bms3.zones, [])

        os.remove('_bms1.json')
        os.remove('_bms1.dat')
        os.remove('_bms1.npy')

class TestCMBlocky(unittest.TestCase):
    def setUp(self):
        cm_geo = mulgrid().rectangular(
            [10.] * 10,
            [10.] * 10,
            [10.] * 10,
            convention=0,
            atmos_type=1, # one per col
            origin=[0,0,0],
            justify='r',
            case=None,
            chars=ascii_lowercase)

        cm_grid = t2grid().fromgeo(cm_geo)
        r = rocktype('ignim')
        cm_grid.add_rocktype(r)
        cm_grid.blocklist[555].rocktype = r

        self.cm_fine = CM_Blocky(cm_geo, cm_grid)

        cm_geo = mulgrid().rectangular(
            [50.] * 2,
            [50.] * 2,
            [50.] * 2,
            convention=0,
            atmos_type=1, # one per col
            origin=[0,0,0],
            justify='r',
            case=None,
            chars=ascii_lowercase)

        cm_grid = t2grid().fromgeo(cm_geo)
        r = rocktype('ignim')
        cm_grid.add_rocktype(r)
        cm_grid.blocklist[6].rocktype = r

        self.cm_coarse = CM_Blocky(cm_geo, cm_grid)

    def test_simple(self):
        bm_geo = mulgrid().rectangular(
            [20.] * 5,
            [20.] * 5,
            [20.] * 5,
            convention=0,
            atmos_type=1, # one per col
            origin=[0,0,0],
            justify='r',
            case=None,
            chars=ascii_lowercase)
        stats, zones = self.cm_fine.populate_model(bm_geo)
        # print(stats, sum(stats))
        stats, zones = self.cm_coarse.populate_model(bm_geo)
        # print(stats, sum(stats))

    def test_shifted(self):
        bm_geo = mulgrid().rectangular(
            [20.] * 5,
            [20.] * 5,
            [20.] * 5,
            convention=0,
            atmos_type=1, # one per col
            origin=[0.1,0.1,0],
            justify='r',
            case=None,
            chars=ascii_lowercase)
        stats, zones = self.cm_fine.populate_model(bm_geo)
        # print(stats, sum(stats))
        stats, zones = self.cm_coarse.populate_model(bm_geo)
        # print(stats, sum(stats), zones)

if __name__ == '__main__':
    unittest.main()
