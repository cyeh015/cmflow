from cmflow.conceptual_models import *
from cmflow.t2data_utils import create_basic_t2data, update_block_geology

import os
import unittest

def test_zonestats_small():
    START = [time.time()]
    geo = mulgrid().rectangular([1000]*20, [1000]*20, [100]*5, origin=[0,0,0],
                                 atmos_type=0)
    stats = BMStats(geo=geo)
    poly = Polygon([
        np.array([1147.7 , 14125.4]),
        np.array([487.0  , 13041.7]),
        np.array([302.0  , 11799.4]),
        np.array([592.7  , 10372.2]),
        np.array([1412.1 , 8389.8]),
        np.array([1993.5 , 7808.4]),
        np.array([2865.8 , 7729.1]),
        np.array([3685.1 , 7544.1]),
        np.array([3975.9 , 6460.4]),
        np.array([4002.3 , 3262.2]),
        np.array([3923.0 , 2020.0]),
        np.array([3685.1 , 1094.9]),
        np.array([5720.3 , 275.5]),
        np.array([6645.4 , 381.3]),
        np.array([7914.1 , 698.4]),
        np.array([9156.3 , 857.0]),
        np.array([9922.8 , 804.2]),
        np.array([11693.7, 1676.4]),
        np.array([12486.6, 2786.5]),
        np.array([12539.5, 4927.4]),
        np.array([13913.9, 7121.2]),
        np.array([13543.9, 8680.6]),
        np.array([11376.5, 8574.9]),
        np.array([9473.5 , 8522.0]),
        np.array([8072.7 , 8944.9]),
        np.array([7121.1 , 10160.7]),
        np.array([6989.0 , 11244.4]),
        np.array([6856.8 , 11905.2]),
        np.array([5641.0 , 13041.7]),
        np.array([4980.2 , 13358.9]),
        np.array([3711.6 , 14204.7]),
        np.array([3262.2 , 14733.3]),
        np.array([1967.1 , 14812.6]),
        np.array([1359.2 , 14389.7]),
        ])
    ztop, zbot = 0.0, -100.0
    cm = CM_Prism('resis', poly, ztop, zbot)
    stats.add_cm(cm)
    ztop, zbot = -200.0, -300.0
    cm = CM_Prism('resis', poly, ztop, zbot)
    stats.add_cm(cm)
    stats.save('tmp_cm.json')
    stats_2 = BMStats(geo=geo)
    stats_2.load('tmp_cm.json', load_geo=False)
    # check if the same
    assert stats.zones == stats_2.zones
    assert np.array_equal(stats.stats, stats_2.stats)
    dat = create_basic_t2data(geo)

    # starting with single rock deflt
    orig_name = dat.grid.rocktypelist[0].name
    dat.grid.rename_rocktype(orig_name, 'NA   ') # not in geological model
    # ATMOS
    for i in range(dat.grid.num_atmosphere_blocks):
        update_block_geology(dat, dat.grid.blocklist[i].name, 'ATMOS')
    # others
    for i in range(dat.grid.num_atmosphere_blocks, dat.grid.num_blocks):
        zi = np.argmax(stats.stats[i,:])
        if stats.stats[i,zi] > 0.0:
            # print CM_0_stats[i,zi], CM_0_stats[i,:]
            new_rname = update_block_geology(dat, dat.grid.blocklist[i].name, stats.zones[0])
    geo.write('gtmp.dat')
    dat.write('tmp.dat')

def test_cm_blocky_full():
    START = [time.time()]

    leapfrog = LeapfrogGM(geometry='gtmp_ay2017_03_6')
    leapfrog.import_leapfrog_csv('grid_gtmp_ay2017_03_6_fit.csv')
    leapfrog.write('_leapfrog.json')
    print_wall_time('loaded leapfrog GM', loop_stop=True)

    geo = mulgrid('gtmp_ay2017_03_6_fit.dat')
    print_wall_time('loaded CM geo', loop_stop=True)

    bm_geo = mulgrid('gtmp_ay2017_05_5a.dat')
    print_wall_time('loaded CM geo and BM geo', loop_stop=True)

    print(('  CM has %i blocks: %s' % (geo.num_blocks, geo.filename)))
    print(('  BM has %i blocks: %s' % (bm_geo.num_blocks, bm_geo.filename)))

    cm = CM_Blocky(geo, leapfrog)
    stats, zones = cm.populate_model(bm_geo)

    print_wall_time('Finished all, total wall time:', total=True)

    print((np.nonzero(stats)))
    np.save('_CM_results.npy', stats)

    with open('_CM_results.json', 'w') as f:
        json.dump({
                    "comments": [
                        "CM (leapfrog): " + leapfrog.import_from,
                        "CM geometry: " + geo.filename,
                        "BM geometry: " + bm_geo.filename,
                    ],
                    "zones": zones,
                    "stats": stats.tolist(),
                  }, f, indent=4, sort_keys=True)


def test_cm_fault_full():
    import glob
    START = [time.time()]

    bm_geo = mulgrid('gtmp_ay2017_05_5a.dat')
    print_wall_time('loaded BM geo', loop_stop=True)
    print(('  BM has %i blocks: %s' % (bm_geo.num_blocks, bm_geo.filename)))

    cm_f = CM_Faults(sorted(glob.glob('*.ts')))
    print_wall_time('loaded faults CM w/ %i faults' % cm_f.num_zones, loop_stop=True)

    stats, zones = cm_f.populate_model(bm_geo)
    print_wall_time('Finished all, total wall time:', total=True)

    save_as = '_CM_faults_results'
    print((np.nonzero(stats)))
    np.save(save_as + '.npy', stats)
    with open(save_as + '.json', 'w') as f:
        json.dump({
                    "comments": [
                        "BM geometry: " + bm_geo.filename,
                    ],
                    "zones": zones,
                    "stats": stats.tolist(),
                  }, f, indent=4, sort_keys=True)


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
