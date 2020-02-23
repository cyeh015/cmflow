from conceptual_models import *

import os
import unittest

class TestBMStats(unittest.TestCase):
    def setUp(self):
        self.geo1 = mulgrid().rectangular(
            [10.] * 3,
            [10.] * 2,
            [10.] * 1,
            convention=0,
            atmos_type=2, # no atm blocks
            origin=[0.,0.,0.],
            justify='r',
            case=None,
            chars=ascii_lowercase)

    def test_empty_save_load(self):
        self.assertEqual(len(self.geo1.block_name_list), 6)

        # new/empty bms requires geo
        with self.assertRaises(BMStatsError):
            bms = BMStats()

        # new/empty bms
        bms = BMStats(geo=self.geo1)
        self.assertEqual(bms.geo, self.geo1)
        self.assertEqual(bms.stats.shape, (6, 0)) # np.array (num_blocks, num_zones)
        self.assertEqual(bms.zones, [])

        # save bms
        self.geo1.write('_bms1.dat')
        bms.save('_bms1.json')
        self.assertTrue(os.path.isfile('_bms1.json'))
        self.assertTrue(os.path.isfile('_bms1.dat'))
        self.assertTrue(os.path.isfile('_bms1.npy'))

        # load bms, with pre-loaded geo
        bms2 = BMStats('_bms1.json', geo=self.geo1)
        self.assertEqual(bms2.geo, self.geo1)
        self.assertEqual(bms2.stats.shape, (6, 0))

        # load bms, will load "geometry" automatically
        bms3 = BMStats('_bms1.json')
        self.assertEqual(bms3.geo.block_name_list, self.geo1.block_name_list)
        self.assertEqual(bms3.stats.shape, (6, 0))
        self.assertEqual(bms3.zones, [])

        os.remove('_bms1.json')
        os.remove('_bms1.dat')
        os.remove('_bms1.npy')

    def test_new(self):
        # new with something with error
        with self.assertRaises(BMStatsError):
            zones = ['A', 'BB', 'CCC']
            stats = np.zeros((6,2))
            bm4 = BMStats(geo=self.geo1, stats=stats, zones=zones)

        # new with something
        zones = ['A', 'BB']
        stats = np.zeros((6,2))
        bm5 = BMStats(geo=self.geo1, stats=stats, zones=zones)
        b1 = self.geo1.block_name_list[0]
        b2 = self.geo1.block_name_list[-1]
        self.assertEqual(list(bm5.cellstats[b1]), [0., 0.])
        self.assertEqual(list(bm5.cellstats[b2]), [0., 0.])

        # add a zone
        bm5.add_zone('CCC', [0., 0., 0., 0., 0., 0.])
        self.assertEqual(list(bm5.cellstats[b1]), [0., 0., 0.])
        self.assertEqual(list(bm5.cellstats[b2]), [0., 0., 0.])

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
        bms_fine = self.cm_fine.calc_bmstats(bm_geo)
        self.assertTrue(np.array_equal(bms_fine.stats, stats))
        self.assertEqual(bms_fine.zones, zones)

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
