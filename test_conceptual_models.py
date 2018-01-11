from conceptual_models import *
import unittest

class TestCM(unittest.TestCase):
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

        self.cm_fine = CM(cm_geo, cm_grid)

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

        self.cm_coarse = CM(cm_geo, cm_grid)

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
        # print stats, sum(stats)
        stats, zones = self.cm_coarse.populate_model(bm_geo)
        # print stats, sum(stats)

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
        print stats, sum(stats)
        stats, zones = self.cm_coarse.populate_model(bm_geo)
        print stats, sum(stats), zones

if __name__ == '__main__':
    unittest.main()
