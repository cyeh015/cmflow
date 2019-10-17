from mulgrids import *
from t2data import *

import numpy as np

# for blocky CM
from shapely.geometry import Polygon
from shapely.geometry import LineString
from rtree import index

# for faults CM
from geom_3dface_utils import Face3D
#
from geom_surface_utils import get_columns_intersect_polygon
from geom_surface_utils import geo_column_polygon

import json
import time

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
        elif grid is None:
            # customised cm
            # !!! NOTE, .populate_model() might not work
            self.zones = []
            self.block = {}
        else:
            raise Exception("Unable to load t2grid or LeapfrogGM object")

        # caching objects, so can call populate multiple times efficiently
        self._inter_areas, self._bm_col_ccis = None, None
        self._inter_lengths, self._bm_lay_clis = None, None
        self._bm_idx = None

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
        if self._inter_areas is not None and self._bm_col_ccis is not None:
            return self._inter_areas, self._bm_col_ccis

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

        self._inter_areas, self._bm_col_ccis = areas, bm_col_ccis
        return areas, bm_col_ccis

    def layer_intersect_length(self, geo):
        if self._inter_lengths is not None and self._bm_lay_clis is not None:
            return self._inter_lengths, self._bm_lay_clis

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

        self._inter_lengths, self._bm_lay_clis = lengths, bm_lay_clis
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
        if self._bm_idx is None:
            self._bm_idx = setup_block_name_index_fast(geo)
        print_wall_time('  created col/lay idx for BM (new method)', loop_stop=True)

        ### calculating and fillinf stats here
        for ii,(bm_ci,bm_li) in enumerate(self._bm_idx):
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
                        zi = self.block[cb]
                        # zi = self.block[self.geo.block_name_index(cb)]
                        stats[ii,zi] = stats[ii,zi] + ivol / bvol
        print_wall_time('  Finished calculating/filling stats', loop_stop=True)
        return stats, self.zones


class CM_Prism(object):
    def __init__(self, name, polygon, ztop, zbottom):
        """ create a conceptual mode
        """
        super(CM_Prism, self).__init__()
        self.name = name
        self.polygon = polygon
        self.ztop, self.zbottom = ztop, zbottom

    def column_intersect_area(self, bm_geo):
        bm_polys = geo_column_polygon(bm_geo)
        areas = np.zeros(bm_geo.num_columns)
        for i,bm_poly in enumerate(bm_polys):
            areas[i] = self.polygon.intersection(bm_poly).area
        print_wall_time('    finished creating column intersection array', loop_stop=True)
        return areas

    def layer_intersect_length(self, bm_geo):
        bm_lines = [LineString([(lay.bottom,0), (lay.top,0)]) for lay in bm_geo.layerlist]
        cm_line = LineString([(self.ztop,0), (self.zbottom,0)])
        lengths = np.zeros(bm_geo.num_layers)
        for i,bm_line in enumerate(bm_lines):
            lengths[i] = bm_line.intersection(cm_line).length
        return lengths

    def populate_model(self, bm_geo):
        """ This is the core of the CM processing, fill-in the stats array.

        stats array rows are of base model blocks, and columns of the zones from
        CM.  Each cell is the portion of the block occupied by the zone.   In
        most cases, the total of each row should be 1.0.
        """
        def setup_block_name_index_fast(bm_geo):
            """ based on mulgrid.setup_block_name_index()

            Note atmosphere blocks may not have proper column or layer index,
            None would be used in place.
            """
            block_ij_list = [] # (coli, layj)
            if bm_geo.num_layers > 0:
                if bm_geo.atmosphere_type  ==  0: # one atmosphere block
                    # bn = bm_geo.block_name(bm_geo.layerlist[0].name, bm_geo.atmosphere_column_name)
                    block_ij_list.append((None, None))
                elif bm_geo.atmosphere_type == 1: # one atmosphere block per column
                    for i,col in enumerate(bm_geo.columnlist):
                        # bn = bm_geo.block_name(bm_geo.layerlist[0].name, col.name)
                        block_ij_list.append((i, None))
                for j,lay in enumerate(bm_geo.layerlist[1:]):
                    for i,col in [(ii,col) for ii,col in enumerate(bm_geo.columnlist) if col.surface > lay.bottom]:
                        # bn = bm_geo.block_name(lay.name, col.name)
                        block_ij_list.append((i, j+1))
            return block_ij_list

        inter_areas = self.column_intersect_area(bm_geo)
        print_wall_time('  column_intersect_area() finished: ', loop_stop=True)
        inter_lengths = self.layer_intersect_length(bm_geo)
        print_wall_time('  layer_intersect_length() finished: ', loop_stop=True)

        stats = np.zeros((bm_geo.num_blocks, 1))

        bm_idx = setup_block_name_index_fast(bm_geo)
        print_wall_time('  created col/lay idx for BM (new method)', loop_stop=True)

        ### calculating and fill-in stats here
        for ii,(bm_ci,bm_li) in enumerate(bm_idx):
            if bm_ci is None or bm_li is None:
                # atmosphere blocks, skip
                continue
            bvol = bm_geo.block_volume(bm_geo.layerlist[bm_li], bm_geo.columnlist[bm_ci])
            ivol = inter_areas[bm_ci] * inter_lengths[bm_li]
            if ivol > 0.0:
                stats[ii,0] = stats[ii,0] + ivol / bvol
        print_wall_time('  Finished calculating/filling stats', loop_stop=True)
        return stats, [self.name]


class CM_Faults(object):
    """ Conceptual Model of faults as 3D surface (Face3D *.ts objects)

    NOTE this simplements the simple way of getting blocks crossed by faults.
    Instead of 3D Face cutting across 3D blocks/elements, I simply let 3D Face
    cuts across layer centre plane.

    if dilation is specified as a positive float, then the line will be dilated
    (.buffer) with the specified amount.  This changes the behaviour of stats,
    instead of two possible values of 0.0/1.0 for normal fault line case, this
    will return stats with intersection area ratio as other CMs.
    """
    def __init__(self, faults=None, dilation=None):
        import os.path
        super(CM_Faults, self).__init__()
        self.zones = []
        self.fault = {}
        self.dilation = None
        if isinstance(dilation, float):
            if dilation > 0.0:
                self.dilation = dilation
                print 'Fault uses dilation'
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

        def column_intersect_area(dilated_line, bm_geo):
            bm_polys = geo_column_polygon(bm_geo)
            areas = np.zeros(bm_geo.num_columns)
            for i,bm_poly in enumerate(bm_polys):
                areas[i] = dilated_line.intersection(bm_poly).area
            print_wall_time('    finished creating column intersection array', loop_stop=True)
            return areas

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
                    block_ij_idx[(None, 0)] = bi
                    bi += 1
                elif geo.atmosphere_type == 1: # one atmosphere block per column
                    for i,col in enumerate(geo.columnlist):
                        # bn = geo.block_name(geo.layerlist[0].name, col.name)
                        block_ij_idx[(i, 0)] = bi
                        bi += 1
                for j,lay in enumerate(geo.layerlist[1:]):
                    for i,col in [(ii,col) for ii,col in enumerate(geo.columnlist) if col.surface > lay.bottom]:
                        # bn = geo.block_name(lay.name, col.name)
                        block_ij_idx[(i, j+1)] = bi
                        bi += 1
            return block_ij_idx

        stats = np.zeros((geo.num_blocks, self.num_zones))
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
                line = LineString([tuple(pt[:2]) for pt in pts])
                if self.dilation is not None:
                    line = line.buffer(self.dilation)
                    for ii in col_idx.intersection(line.bounds):
                        iarea = line.intersection(col_polys[ii]).area
                        if iarea > 0.0:
                            if (ii,jj) in block_ij_idx:
                                bi = block_ij_idx[(ii,jj)]
                                c = geo.column_name(geo.block_name_list[bi])
                                carea = geo.column[c].area
                                stats[bi,fi] = stats[bi,fi] + iarea / carea
                                count += 1
                else:
                    for ii in col_idx.intersection(line.bounds):
                        if line.intersection(col_polys[ii]).length > 0.0:
                            if (ii,jj) in block_ij_idx:
                                bi = block_ij_idx[(ii,jj)]
                                stats[bi,fi] = 1.0
                                count += 1
                                # print 'found ', geo.block_name_list[bi]
            print_wall_time('  finished layer %i, found %i blocks' % (jj, count), loop_stop=True)
        return stats, self.zones


class BMStats(object):
    """ Base Model Stats, mainly numpy arrays with rows corresponding to mulgrid
    blocks, and columns corresponding to zones.  Each is a value, usually
    between 0.0 and 1.0.  Often 1.0 is indicating that particular block is fully
    within the zone.

    .stats numpy array (n,m), n = num of model blocks, m = num of zones
    .zones list of zone names (str)
    .zonestats dict of stats column by zone names
    .cellstats dict of stats row by block name
    """
    def __init__(self, geo):
        self.bmgeo = geo
        n = geo.num_blocks
        self.stats = np.zeros((n, 0)) # 'empty' array, ready to concatenate etc
        self.zones = []
        self._reindex()

    def _reindex(self):
        self.zonestats = {z:self.stats[:,i] for i,z in enumerate(self.zones)}
        self.cellstats = {b:self.stats[i,:] for i,b in enumerate(self.bmgeo.block_name_list)}

    def save(self, filename):
        import os.path
        root, ext = os.path.splitext(filename)
        npy_filename = root + '.npy'
        np.save(npy_filename, self.stats)
        # the .npy file will be expected in the same directory, so strip dir
        npy_filename = os.path.split(npy_filename)[1]
        with open(filename, 'w') as f:
            json.dump({
                        "comments": [
                            "BM geometry: %s" % self.bmgeo.filename,
                        ],
                        "stats": npy_filename,
                        "zones": self.zones,
                      }, f, indent=True, sort_keys=True)

    def load(self, filename):
        import os.path
        with open(filename, 'r') as f:
            data = json.load(f)
            self.zones = data['zones']
            # .npy file expected relative to the json file
            npy_filename = os.path.join(os.path.split(filename)[0], data["stats"])
            self.stats = np.load(npy_filename)
            n = self.stats.shape[0]
            if n != self.bmgeo.num_blocks:
                msg1 = 'Loaded BMStats has different number of blocks to geometry file.'
                msg2 = 'BMStats (%i) != Geometry (%i)' % (n, self.bmgeo.num_blocks)
                raise Exception('\n'.join([msg1, msg2]))
            self._reindex()

    def add_stats(self, stats, zones):
        for i,zz in enumerate(zones):
            ss = stats[:,i:i+1]
            if zz in self.zones:
                ii = self.zones.index(zz)
                self.stats[:,ii] = self.stats[:,ii] + ss[:,0]
            else:
                self.stats = np.concatenate((self.stats, ss), axis=1)
                self.zones.append(zz)
        self._reindex()

    def add_cm(self, cm):
        stats, zones = cm.populate_model(self.bmgeo)
        self.add_stats(stats, zones)

def test_zonestats_small():
    from t2data_utils import create_basic_t2data, update_block_geology
    START = [time.time()]
    geo = mulgrid().rectangular([1000]*20, [1000]*20, [100]*5, origin=[0,0,0],
                                 atmos_type=0)
    stats = BMStats(geo)
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
    stats_2 = BMStats(geo)
    stats_2.load('tmp_cm.json')
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
    test_zonestats_small()
    pass