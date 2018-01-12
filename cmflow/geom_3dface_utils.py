
def line_intersect_plane(p0, p1, pp, pn, in_segment=True,
    check_on_plane=True, tol=0.000001):
    """ Return a point (tuple of three floats) in space where line segment p0-p1
    intersects the plane defined by point-normal form (pp, pn).

    Note, set both options to True is most desirable if what to know exactly how
    the segment intersect the plane.
    --- complicated cases:

    1. Return None if 'check_on_plane' is False AND if segment is on the plane.

    2. Return (0, 1) if 'check_on_plane' is True AND if segment is on the
    plane.

    3. Return 0 or 1 if either p0 or p1 is on the plane when 'in_segment' is
    True.

    4. Return None if segment is parallel to the plane but not on the plane.

    5. Return None if segment's extension intersects the plane, but not between
    the two points when 'in_segment' is True.

    6. Return a valid point if segment or its extension intersects the plane
    when 'in_segment' is False.
    """
    def add_3(a, b):
        return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    def sub_3(a, b):
        return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    def mul_3(a, f):
        return tuple([x*f for x in a])
    def cross_3(a, b):
        """ cross product of two vectors in 3D space """
        return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    def dot_3(a, b):
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

    #####
    u = sub_3(p1, p0)
    pn_dot_line = dot_3(pn, u) # dot prod of line and plane normal

    if abs(pn_dot_line) > tol:
        # the factor of the point between p0 -> p1 (0 - 1)
        w = sub_3(p0, pp)
        fac = -dot_3(pn, w) / pn_dot_line
        # if 'fac' is between (0 - 1) the point intersects with the segment.
        # otherwise:
        #  < 0.0: behind p0.
        #  > 1.0: infront of p1.
        if in_segment:
            if abs(fac - 0.0) < tol:
                # intersect at p0
                return 0
            elif abs(fac - 1.0) < tol:
                # intersect at p1
                return 1
            elif fac < 0.0 or fac > 1.0:
                return None
        pi = add_3(p0, mul_3(u, fac))
        return pi
    else:
        # The segment is parallel to plane
        if check_on_plane:
            v = sub_3(p1, pp)
            pn_dot_pp1 = dot_3(pn, v)
            if abs(pn_dot_pp1) < tol:
                # the segment is ON the plane
                return (0, 1)
        return None


def point_on_plane(p, pp, pn, tol=0.00001):
    def add_3(a, b):
        return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    def sub_3(a, b):
        return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    def mul_3(a, f):
        return tuple([x*f for x in a])
    def cross_3(a, b):
        """ cross product of two vectors in 3D space """
        return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])
    def dot_3(a, b):
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

    #####
    u = sub_3(p, pp)
    pn_dot_line = dot_3(pn, u) # dot prod of line and plane normal

    if abs(pn_dot_line) < tol:
        return True
    else:
        return False

class Point(object):
    def __init__(self, name='', coordinates=[]):
        self.name = name
        self.coordinates = coordinates

        # some back references for searching etc
        self.in_segments = []
        self.in_triangles = []

    @property
    def co(self):
        """ short hand for coordinates"""
        return self.coordinates

    def __repr__(self):
        return 'Point:'+self.name

class Triangle(object):
    def __init__(self, points=[]):
        self.points = points
        self.index = None

        # supplemental data
        self.segments = []
        self.intersections = []

    def edges(self, excludes=[]):
        """ excludes is a list of points/segments.  Points are single values,
        and segments are tuples. """
        from itertools import combinations
        all_verts = [x for x in self.points if x not in excludes]
        all_combs = []
        for p0,p1 in combinations(all_verts, 2):
            if (p0,p1) not in excludes and (p1,p0) not in excludes:
                all_combs.append((p0,p1))
        return all_combs

    def __repr__(self):
        return 'Triangle:' + ' '.join([p.name for p in self.points])

class Face3D(object):
    def __init__(self):
        self.header = {}
        self.points = {}
        self.triangles = []

        # utility storage for faster searching etc
        # self.setup_supplemental_mappings()

    # def setup_supplemental_mappings(self):
    #     self.points_in_triangles = {}

    #     self.triangle_segments = []
    #     self._seg_intersections = [] # each is two list ([cross coordinates],[cross corner name])
    #     self.seg_intersections = {} # keyed with two-way tuples, (a,b) and (b,a) point to the same object in list

    #     for tri in self.triangles:
    #         from itertools import combinations
    #         for p0,p1 in combinations(tri, 2):
    #             self._seg_intersections.append(None)
    #             seg_i = len(self._seg_intersections)
    #             self.seg_intersections[(p0,p1)] = seg_i
    #             self.seg_intersections[(p1,p0)] = seg_i
    #             self.triangle_segments.append(seg_i)

    @property
    def num_triangles(self):
        return len(self.triangles)

    def read(self, filename):
        # TODO: implement full support, see http://paulbourke.net/dataformats/gocad/gocad.pdf
        f = open(filename, 'r')
        if not f.next().startswith('GOCAD TSurf 1'):
            raise Exception("File type not supported, need .ts file starts with 'GOCAD TSurf 1'.")
        done = False
        for line in f:
            if line.startswith('HEADER'):
                line = f.next()
                while '}' not in line:
                    c = line.strip().split(':')
                    self.header[c[0]] = c[1]
                    line = f.next()
            elif line.startswith('TFACE'):
                print 'Reading TFACE...',
                for line in f:
                    if line.startswith('PVRTX'):
                        c = line.strip().split()
                        self.points[c[1]] = Point(c[1], [float(x) for x in c[2:]])
                    if line.startswith('TRGL'):
                        c = line.strip().split()
                        self.triangles.append(Triangle([self.points[p] for p in c[1:]]))
                        # add back references
                        self.triangles[-1].index = len(self.triangles) - 1
                        for p in c[1:]:
                            self.points[p].in_triangles.append(self.triangles[-1])
                    if line.startswith('END'):
                        done = True
                        break
            else:
                pass
        f.close()
        if done:
            print 'Done.'

    def __repr__(self):
        if 'name' in self.header:
            name = '%s - ' % self.header['name']
        return name + 'TSurf with %i points and %i triangles' % (len(self.points),len(self.triangles))

    def set_cutting_plane(self, pl_point, pl_normal):
        # cached result of line-plane intersections, (a,b) and (b,a) keeps the same result
        self.cached_seg_intersections = {} # key is (Point,Point)
        self.cached_tri_intersections = [] # list same order as self.triangles

        self._calc_all_seg_intersections(pl_point, pl_normal)
        self._proc_all_tri_intersections(pl_point, pl_normal)

        # print '+ number of tri(0,2):', len(self.find_tri_match_itype((0,2)))
        # print '+ number of tri(1,0):', len(self.find_tri_match_itype((1,0)))
        # print '+ number of tri(2,0):', len(self.find_tri_match_itype((2,0)))
        # print '+ number of tri(3,0):', len(self.find_tri_match_itype((3,0)))
        # print '+ number of tri(1,1):', len(self.find_tri_match_itype((1,1)))
        # print '+ number of tri(0,0):', len(self.find_tri_match_itype((0,0)))
        # print '- total two pointers:', len(self.find_tri_match_itype((0,2))) + len(self.find_tri_match_itype((2,0))) + len(self.find_tri_match_itype((1,1)))

    def search_line(self):
        """ this traverses along triangles that have total of two intersections
        with plane, it only goes through triangle that has two intersection
        points (2,0), (0,2), or (1,1) """

        # these are all storing actual ref to triangle objects
        self.cached_tri_traversed = []
        self.cached_tri_on_plane = self.find_tri_match_itype((3,0))
        self.cached_tri_one_corner = self.find_tri_match_itype((1,0))

        # starting triangle must be valid (two intersections)
        cur_tri = None
        for tri_i, isect in enumerate(self.cached_tri_intersections):
            if sum(isect['itype']) == 2 :
                # tri_start will be usable starting point, (3,0) or (0,0) is not useful, (1,0) also not ideal
                cur_tri = self.triangles[tri_i]
                break

        if cur_tri is None:
            return []

        self.cached_tri_traversed.append(cur_tri)

        isect = self.cached_tri_intersections[cur_tri.index]
        # cur_ep can be either a segment or a point
        if isect['itype'] == (0,2):
            # two edges crosses plane
            ks = isect['crosses'].keys()
            cur_ep, other_end_ep = ks[0], ks[1]
        elif isect['itype'] == (2,0):
            # two corners on plane
            cur_ep, other_end_ep = isect['corners'][0], isect['corners'][1]
        else:
            # one corner on plane, one edge crosses plane
            ks = isect['crosses'].keys()
            cur_ep, other_end_ep = isect['corners'][0], ks[0]

        line = [self.get_intersection_co(other_end_ep, isect), self.get_intersection_co(cur_ep, isect)]
        line += self._traverse_one_way(cur_tri, cur_ep, other_end_ep)
        line = line[::-1]
        line += self._traverse_one_way(cur_tri, other_end_ep, other_end_ep)

        # print 'finished, len=', len(line)
        return line
        # print '\n'.join([('%.2f, %.2f, %.2f' % tuple(coo)) for coo in line])

    def _traverse_one_way(self, cur_tri, cur_ep, cur_ep_corner2):
        """ Note that 'cur_ep_corner2' is for special cases where the starting
        triangle has intersection of type (2,0) and finding neighbours using
        only one point, this gets rid of the neighbour that shares exactly the
        same edge that is on the plane. """
        thru_coordinates = []
        while True:
            if isinstance(cur_ep, tuple):
                # sharing edge
                nbs = self._neighbours_by_edge(cur_tri, cur_ep,
                    self.cached_tri_traversed + self.cached_tri_on_plane + self.cached_tri_one_corner)
                # if there is good neighbours at all, it must be useful

            else:
                # corner
                nbs = self._neighbours_by_corner(cur_tri, cur_ep,
                    self.cached_tri_traversed + self.cached_tri_on_plane + self.cached_tri_one_corner)
                # can have good neighbour that share same cur_ep,
                nbs_ = [nb for nb in nbs]
                for nb in nbs_:
                    isect = self.cached_tri_intersections[nb.index]
                    if isect['itype'] == (2,0):
                        if cur_ep in isect['corners'] and cur_ep_corner2 in isect['corners']:
                            nbs.remove(nb)

            if len(nbs) == 0:
                break

            # if len(nbs) != 1:
            #     print '!!!',  cur_ep, cur_tri, 'has nbs:', [x for x in nbs], [self.cached_tri_intersections[x.index]['itype'] for x in nbs]
            #     raise Exception('Triangle %i has %i neighbours.' % (cur_tri.index, len(nbs)))

            isect = self.cached_tri_intersections[nbs[0].index]
            exits = []
            for ex in isect['crosses'].keys() + isect['corners']:
                if isinstance(cur_ep, tuple):
                    if ex in [cur_ep, cur_ep[::-1]]:
                        continue
                else:
                    if ex == cur_ep:
                        continue
                exits.append(ex)
            # exits = list(set(all_eps) - set([cur_ep]))
            if isect['itype'] == (2,0) and not isinstance(cur_ep, tuple):
                # only if traversing through a (2,0) triangle (ie. from one corner to another)
                cur_ep_corner2 = cur_ep
            else:
                cur_ep_corner2 = None

            if len(exits) == 0:
                break

            # if len(exits) != 1:
            #     print '$$$', nbs[0], exits, cur_ep, self.cached_tri_intersections[nbs[0].index]['itype']
            #     raise Exception('Triangle %i has %i exits.' % (nbs[0].index, len(exits)))

            # update current
            cur_tri, cur_ep = nbs[0], exits[0]
            self.cached_tri_traversed.append(cur_tri)
            thru_coordinates.append(self.get_intersection_co(cur_ep, isect))
            # print cur_tri.index, self.cached_tri_intersections[cur_tri.index]['itype'] # , thru_coordinates[-1]

        return thru_coordinates

    def _neighbours_by_edge(self, tri, seg, excludes=[]):
        if isinstance(tri, int):
            tri = self.triangles[tri]
        excludes.append(tri)
        tri_sets = [set(p.in_triangles) for p in seg]
        nbs = tri_sets[0]
        for s in tri_sets[1:]:
            nbs.intersection_update(s)
        return list(nbs - set(excludes))

    def _neighbours_by_corner(self, tri, p, excludes=[]):
        if isinstance(tri, int):
            tri = self.triangles[tri]
        excludes.append(tri)
        nbs = []
        for t in p.in_triangles:
            if t not in excludes:
                nbs.append(t)
        return nbs

    def _calc_all_seg_intersections(self, pl_point, pl_normal):
        for i, tri in enumerate(self.triangles):
            for p0, p1 in tri.edges():
                if (p0,p1) not in self.cached_seg_intersections:
                    isect = line_intersect_plane(p0.co, p1.co, pl_point, pl_normal,
                        in_segment=True, check_on_plane=True)

                    # translate to Point references if needed
                    if isect == 0:
                        isect = [p0]
                    elif isect == 1:
                        isect = [p1]
                    elif isect == (0,1):
                        isect = [p0,p1]
                    else:
                        # just leave it as coordinates [x,y,z] or None
                        pass
                    self.cached_seg_intersections[(p0,p1)] = isect
                    self.cached_seg_intersections[(p1,p0)] = isect
        # NOTE: isect has a few possible values:
        # None - no intersection
        # 0 - p0 on plane  ->  tranlated to p0
        # 1 - p1 on plane  ->  tranlated to p1
        # (0,1) - line on plane  ->  tranlated to (p0,p1)
        # [x,y,z] - intersect in the segment

    def get_intersection_co(self, ep, isect):
        """ ep is either a Point or an edge (tuple of Point), isect is the
        intersection info dict """
        if isinstance(ep, tuple):
            return isect['crosses'][ep]
        else:
            return ep.co

    def find_tri_match_itype(self, itype):
        """ create a list of triangles that qualify the itype """
        return [self.triangles[i] for i,v in enumerate(self.cached_tri_intersections) if v['itype'] == itype]

    def _proc_all_tri_intersections(self, pl_point, pl_normal):
        # go thru all triangles, and determine which case it is
        for tri in self.triangles:
            corners = []
            crosses = {}
            for seg in tri.edges():
                isect = self.cached_seg_intersections[seg]
                if isect is None:
                    continue
                elif len(isect) <= 2:
                    corners += isect
                else:
                    # three values, must be [x,y,z]
                    crosses[seg] = isect

            corners = list(set(corners))
            self.cached_tri_intersections.append({
                'itype': (len(corners),len(crosses)),
                'corners': corners,
                'crosses': crosses,
                })

