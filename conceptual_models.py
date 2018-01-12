from mulgrids import *
from t2data import *

# for blocky CM
from shapely.geometry import Polygon
from shapely.geometry import LineString
from rtree import index

# for faults CM
from geom_3dface_utils import Face3D

import json
import time

import matplotlib.pyplot as plt
from wairakei_res import draw_res, insert_res_as_wells

START = [time.time()]
def print_wall_time(msg, loop_stop=False, total=False):
    t = time.time()
    if total:
        print msg, ': %f seconds' % (t - START[0])
    else:
        print msg, ': %f seconds' % (t - START[-1])
    if loop_stop:
        START.append(t)

class LeapfrogGM(object):
    def __init__(self, geometry=''):
        super(LeapfrogGM, self).__init__()
        self.litholist = []
        self.blocklitho = {}
        self.import_from = '' # optional, if imported from leapfrog
        self.geometry = '' # optional, matching mulgrid geometry file

    def import_leapfrog_csv(self, filename, report=False):
        """ load geology info from Leapfrog's 'Generate rock types' feature.  The
        .csv file from Leapfrog usually starts with a table of 'LithoCode,Lithology'
        then a longer table of 'BlockName,LithoCode'.  Both of these will be
        retuened as dictionaries. """
        import csv
        f = open(filename,'r')
        allrows = csv.reader(f)
        lithocodes, self.litholist = {}, []
        self.blocklitho = {}
        # TODO: this file reading is ugly, need work
        read_litho_name, read_block_litho = False, False
        for row in allrows:
            if len(row) == 1: continue
            if len(row) == 3:
                if row[0] == '#' and row[1] == 'LithoCode':
                    read_litho_name = True
                    continue
            if len(row) == 2:
                if row[0] == 'BlockName' and row[1] == 'LithoCode':
                    read_litho_name = False
                    read_block_litho = True
                    continue
            if read_litho_name:
                if not row[0] == '#': raise Exception
                self.litholist.append(row[2])
                lithocodes[int(row[1])] = len(self.litholist)-1
            if read_block_litho:
                if len(row) <> 2: break
                self.blocklitho[row[0]] = lithocodes[int(row[1])]
        f.close()
        self.import_from = filename
        if report:
            print '%8s Lithology found.' % len(self.litholist)
            print '%8s Blocks allocated.' % len(self.blocklitho)

    def write(self, filename):
        with open(filename, 'w') as f:
            data = {
                "comments": [
                    "Simplified Leapfrog Geological Model",
                    "zones: a list of zone names",
                ],
                "import_from" : self.import_from,
                "geometry" : self.geometry,
                "zones": self.litholist,
                "blocks": self.blocklitho,
            }
            json.dump(data, f, indent=4, sort_keys=True)

    def read(self, filename):
        with open(filename, 'r') as f:
            data = json.load(f)
            self.import_from = data['import_from']
            self.geometry = data['geometry']
            self.litholist = data['zones']
            self.blocklitho = data['blocks']

class CM_Blocky(object):
    """ A nceptual Model object is a model that represents a space by a list of
    zones (usually exists as a blocky model, each block has a *rocktype*). """
    def __init__(self, geo, grid):
        """ initialise CM model by load mulgrid geometry and t2grid objects """
        super(CM_Blocky, self).__init__()

        if isinstance(geo, str):
            self.geo = mulgrid(geo)
        elif isinstance(geo, mulgrid):
            self.geo = geo
        else:
            raise Exception("Unable to load mulgrid geometry file")

        if isinstance(grid, str):
            self._load_from_t2grid(t2data(grid).grid)
        elif isinstance(grid, t2data):
            self._load_from_t2grid(grid.grid)
        elif isinstance(grid, t2grid):
            self._load_from_t2grid(grid)
        elif isinstance(grid, LeapfrogGM):
            # if LeapfrogGM object
            self.zones = grid.litholist
            self.block = grid.blocklitho
        else:
            raise Exception("Unable to load t2grid or LeapfrogGM object")

        self.num_zones = len(self.zones)

    def _load_from_t2grid(self, grid):
        print_wall_time('_load_from_t2grid()', loop_stop=True)
        self.zones = [r.name for r in grid.rocktypelist]
        print_wall_time('  created zones from t2grid', loop_stop=True)
        self.block = {b: grid.rocktypelist.index(grid.block[b].rocktype) for b in self.geo.block_name_list}
        print_wall_time('  created block zone dict (old method)', loop_stop=True)

    def column_intersect_area(self, geo):
        """ return an two-d array of (absolute) area of column intersections
        between the model column and CM columns.

        'bm_col_ccis' is a list the same length/order as geo.block_name_list, each
        element is a (varying length) list of cm columns that intersects bm
        column.
        """

        bm_polys = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]
        print_wall_time('    constructed all %i BM polys' % geo.num_columns, loop_stop=True)

        # CM usually has larger number of columns and is regular, so I should
        # probably do RTree on CM grid
        cm_idx = index.Index()
        cm_polys = [Polygon([n.pos for n in c.node]) for c in self.geo.columnlist]
        print_wall_time('    constructed all %i CM polys' % self.geo.num_columns, loop_stop=True)
        for i,poly in enumerate(cm_polys):
            cm_idx.insert(i, poly.bounds)
        print_wall_time('    constructed CM polys RTree', loop_stop=True)

        bm_col_ccis = [] # list of CM col indices that intersects BM columns
        areas = np.zeros((geo.num_columns, self.geo.num_columns))
        for i,bm_poly in enumerate(bm_polys):
            bm_col_ccis.append([])
            for j in cm_idx.intersection(bm_poly.bounds):
                areas[i,j] = bm_poly.intersection(cm_polys[j]).area
                if areas[i,j] > 0.0:
                    bm_col_ccis[-1].append(j)
        print_wall_time('    finished creating column intersection array', loop_stop=True)
        return areas, bm_col_ccis

    def layer_intersect_length(self, geo):
        bm_lines = [LineString([(lay.bottom,0), (lay.top,0)]) for lay in geo.layerlist]
        cm_lines = [LineString([(lay.bottom,0), (lay.top,0)]) for lay in self.geo.layerlist]
        bm_lay_clis = [] # list of CM lay indices that intersects BM layers
        lengths = np.zeros((geo.num_layers, self.geo.num_layers))
        for i,bm_line in enumerate(bm_lines):
            bm_lay_clis.append([])
            for j,cm_line in enumerate(cm_lines):
                lengths[i,j] = bm_line.intersection(cm_line).length
                if lengths[i,j] > 0.0:
                    bm_lay_clis[-1].append(j)
        return lengths, bm_lay_clis

    def populate_model(self, geo):
        """ This is the core of the CM processing, fill-in the stats array.

        stats array rows are of base model blocks, and columns of the zones from
        CM.  Each cell is the portion of the block occupied by the zone.   In
        most cases, the total of each row should be 1.0.
        """
        inter_areas, bm_col_ccis = self.column_intersect_area(geo)
        print_wall_time('  column_intersect_area() finished: ',
                        loop_stop=True)
        inter_lengths, bm_lay_clis = self.layer_intersect_length(geo)
        print_wall_time('  layer_intersect_length() finished: ',
                        loop_stop=True)

        stats = np.zeros((geo.num_blocks, self.num_zones))
        print_wall_time('  created stats array %i x %i' % (geo.num_blocks, self.num_zones), loop_stop=True)

        def setup_block_name_index_fast(geo):
            """ based on mulgrid.setup_block_name_index()

            Note atmosphere blocks may not have proper column or layer index,
            None would be used in place.
            """
            block_ij_list = [] # (coli, layj)
            if geo.num_layers > 0:
                if geo.atmosphere_type  ==  0: # one atmosphere block
                    # bn = geo.block_name(geo.layerlist[0].name, geo.atmosphere_column_name)
                    block_ij_list.append((None, None))
                elif geo.atmosphere_type == 1: # one atmosphere block per column
                    for i,col in enumerate(geo.columnlist):
                        # bn = geo.block_name(geo.layerlist[0].name, col.name)
                        block_ij_list.append((i, None))
                for j,lay in enumerate(geo.layerlist[1:]):
                    for i,col in [(ii,col) for ii,col in enumerate(geo.columnlist) if col.surface > lay.bottom]:
                        # bn = geo.block_name(lay.name, col.name)
                        block_ij_list.append((i, j+1))
            return block_ij_list

        # cm_idx = setup_block_name_index_fast(self.geo)
        # print_wall_time('created col/lay idx for CM (new method)', loop_stop=True)
        bm_idx = setup_block_name_index_fast(geo)
        print_wall_time('  created col/lay idx for BM (new method)', loop_stop=True)

        ### calculating and fillinf stats here
        for ii,(bm_ci,bm_li) in enumerate(bm_idx):
            if bm_ci is None or bm_li is None:
                # atmosphere blocks, skip
                continue
            bvol = geo.block_volume(geo.layerlist[bm_li], geo.columnlist[bm_ci])
            # only do check columns and layers that actually intersect current BM block
            for cm_ci in bm_col_ccis[bm_ci]:
                for cm_li in bm_lay_clis[bm_li]:
                    ivol = inter_areas[bm_ci,cm_ci] * inter_lengths[bm_li,cm_li]
                    if ivol > 0.0:
                        cb = self.geo.block_name(self.geo.layerlist[cm_li].name,
                                                 self.geo.columnlist[cm_ci].name)
                        if cb not in self.block:
                            continue
                        # print 'found something CM block %i within BM' % self.geo.block_name_index[cb]
                        # print inter_areas[bm_ci,cm_ci], inter_lengths[bm_li,cm_li]
                        # print bm_ci, bm_col_ccis[bm_ci]
                        # print bm_li, bm_lay_clis[bm_li]

                        # plt.clf()
                        # geo.layer_plot(plt=plt, linecolour='b', column_names=[geo.columnlist[bm_ci].name])
                        # self.geo.layer_plot(plt=plt, linecolour='r', column_names=[self.geo.columnlist[cm_ci].name])
                        # draw_res(plt)
                        # plt.show()

                        zi = self.block[cb]
                        # zi = self.block[self.geo.block_name_index(cb)]
                        stats[ii,zi] = stats[ii,zi] + ivol / bvol
        print_wall_time('  Finished calculating/filling stats', loop_stop=True)
        return stats, self.zones

class CM_Faults(object):
    """ Conceptual Model of faults as 3D surface (Face3D *.ts objects)

    NOTE this simplements the simple way of getting blocks crossed by faults.
    Instead of 3D Face cutting across 3D blocks/elements, I simply let 3D Face
    cuts across layer centre plane.
    """
    def __init__(self, faults=None):
        import os.path
        super(CM_Faults, self).__init__()
        self.zones = []
        self.fault = {}
        if faults is None:
            pass
        elif isinstance(faults, list):
            # a list of *.ts files to load
            for filename in faults:
                fault = Face3D()
                fault.read(filename)
                zonename = os.path.splitext(filename)[0]
                self.zones.append(zonename)
                self.fault[zonename] = fault
        elif isinstance(faults, dict):
            # dictionary of Face3D objects
            for zonename in sorted(faults.keys()):
                self.zones.append(zonename)
                self.fault = faults
        else:
            raise Exception
        self.num_zones = len(self.zones)

    def populate_model(self, geo):

        def column_polygons(geo):
            # CM usually has larger number of columns and is regular, so I should
            # probably do RTree on CM grid
            column_idx = index.Index()
            column_polys = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]
            print_wall_time('    constructed all %i column polygons' % geo.num_columns, loop_stop=True)
            for i,poly in enumerate(column_polys):
                column_idx.insert(i, poly.bounds)
            print_wall_time('    constructed column polygons RTree', loop_stop=True)
            return column_polys, column_idx

        def setup_block_name_index_fast(geo):
            """ based on mulgrid.setup_block_name_index()

            Note atmosphere blocks may not have proper column or layer index,
            None would be used in place.
            """
            block_ij_idx = {} # { (coli, layj): block index }
            bi = 0
            if geo.num_layers > 0:
                if geo.atmosphere_type  ==  0: # one atmosphere block
                    # bn = geo.block_name(geo.layerlist[0].name, geo.atmosphere_column_name)
                    block_ij_idx[(None, None)] = bi
                    bi += 1
                elif geo.atmosphere_type == 1: # one atmosphere block per column
                    for i,col in enumerate(geo.columnlist):
                        # bn = geo.block_name(geo.layerlist[0].name, col.name)
                        block_ij_idx[(i, None)] = bi
                        bi += 1
                for j,lay in enumerate(geo.layerlist[1:]):
                    for i,col in [(ii,col) for ii,col in enumerate(geo.columnlist) if col.surface > lay.bottom]:
                        # bn = geo.block_name(lay.name, col.name)
                        block_ij_idx[(i, j+1)] = bi
                        bi += 1
            return block_ij_idx

        stats = np.zeros((geo.num_blocks, self.num_zones), dtype=bool)
        col_polys, col_idx = column_polygons(geo)
        block_ij_idx = setup_block_name_index_fast(geo)
        print_wall_time('  setup_block_name_index_fast()', loop_stop=True)
        for jj,lay in enumerate(geo.layerlist):
            count = 0 # intersected blocks count, per layer
            z = lay.centre
            for fi,fname in enumerate(self.zones):
                fault = self.fault[fname]
                fault.set_cutting_plane((0.0,0.0,z), (0.0,0.0,1.0))
                pts = fault.search_line()
                if len(pts) < 2:
                    # not enough points to construct LineString, layer not
                    # cutting through the fault Face
                    print "    skipping layer %i '%s' with fault '%s'" % (jj, lay.name, fname)
                    continue
                line = LineString([tuple(pt[:2]) for pt in fault.search_line()])
                for ii in col_idx.intersection(line.bounds):
                    if line.intersection(col_polys[ii]).length > 0.0:
                        if (ii,jj) in block_ij_idx:
                            bi = block_ij_idx[(ii,jj)]
                            stats[bi,fi] = True
                            count += 1
                            # print 'found ', geo.block_name_list[bi]
            print_wall_time('  finished layer %i, found %i blocks' % (jj, count), loop_stop=True)
        return stats, self.zones

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

    print '  CM has %i blocks: %s' % (geo.num_blocks, geo.filename)
    print '  BM has %i blocks: %s' % (bm_geo.num_blocks, bm_geo.filename)

    cm = CM_Blocky(geo, leapfrog)
    stats, zones = cm.populate_model(bm_geo)

    print_wall_time('Finished all, total wall time:', total=True)

    print np.nonzero(stats)
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
    print '  BM has %i blocks: %s' % (bm_geo.num_blocks, bm_geo.filename)

    cm_f = CM_Faults(sorted(glob.glob('*.ts')))
    print_wall_time('loaded faults CM w/ %i faults' % cm_f.num_zones, loop_stop=True)

    stats, zones = cm_f.populate_model(bm_geo)
    print_wall_time('Finished all, total wall time:', total=True)

    save_as = '_CM_faults_results'
    print np.nonzero(stats)
    np.save(save_as + '.npy', stats)
    with open(save_as + '.json', 'w') as f:
        json.dump({
                    "comments": [
                        "BM geometry: " + bm_geo.filename,
                    ],
                    "zones": zones,
                    "stats": stats.tolist(),
                  }, f, indent=4, sort_keys=True)


if __name__ == '__main__':
    # test_cm_blocky_full()
    # test_cm_fault_full()
    pass