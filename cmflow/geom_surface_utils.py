from mulgrids import *
from matplotlib import pyplot as plt
from descartes import PolygonPatch

from shapely.geometry import shape
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry import Polygon

import json
import unittest
import os.path

BLUE = '#6699cc'
GRAY = '#999999'
BLACK = '#000000'

def load_feature(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
        if data['type'] == 'Feature':
            g = shape(data['geometry'])
            meta = {
                'filename': filename,
                'type': data['type'],
                'properties': data['properties'],
                'id': data['id']
            }
            return g, meta

def geo_column_polygon(geo, col_name='', cache_shapely=True):
    """ Returns a shapely Polygon object of the column specified by name.

    If col_name is specified (str) then a single polygon will be returned,
    otherwise a list of (all columns') polygons will be returned. If
    cache_shapely is True, it will attempt to use the cached Polygon within geo.

    If caching is on, geo object will have two additional properties, analog to
    .columnlist and .column:
        ._columnpolygonlist
        ._columnpolygon
    """
    if cache_shapely:
        if hasattr(geo, '_columnpolygonlist'):
            c_polygons = geo._columnpolygonlist
        else:
            c_polygons = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]
            geo._columnpolygonlist = c_polygons
            geo._columnpolygon = {c.name: p for c,p in zip(geo.columnlist, geo._columnpolygonlist)}
    else:
        c_polygons = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]
    if col_name:
        return c_polygons[geo.columnlist.index(geo.column[col_name])]
    else:
        return c_polygons

def get_columns_intersect_polygon(polygon, geo,
                                  threshold=0.0,
                                  cache_shapely=True):
    """ finds all columns that intersects shapely Polygon objects of each
    column.  Threshold can be set to ensure intersection is of significant
    portion (of column area that is inside the polygon, 1.0 means fully enclosed
    by the polygon object).

    NOTE column Polygon objects can be cached, so we don't need to create
    shapely Polygon objects of all columns every time.

    Returns a tuple of:
        - column names of columns (that intersects)
        - area portion of column (0.0-1.0) (that intersects)
        - shapely polygons of these columns
    """
    c_polygons = geo_column_polygon(geo, cache_shapely=cache_shapely)

    inter_colnames, inter_portions, column_polygons = [], [], []
    for c,cp in zip(geo.columnlist, c_polygons):
        if polygon.intersects(cp):
            portion = polygon.intersection(cp).area / cp.area
            if portion > threshold:
                inter_colnames.append(c.name)
                inter_portions.append(portion)
                column_polygons.append(cp)
    return inter_colnames, inter_portions, column_polygons

def snap(sorted_vals, value, direction=None, minimum=None, allow_outside=False):
    """ snap a value into one of the provided values (assume sorted low to high)

    Sanpping direction can be None, 'down', or 'up'.  If minimum limit is
    supplied, it will then be moved up until condition is satisfied.
    """
    if direction not in [None, 'down', 'up']:
        raise Exception("Direction of snapping has to be None, 'down', or 'up'")
    if allow_outside:
        if value < sorted_vals[0]:
            print "value %f outside range, snapped to %f" % (value, sorted_vals[0])
            return sorted_vals[0]
        if value > sorted_vals[-1]:
            print "value %f outside range, snapped to %f" % (value, sorted_vals[-1])
            return sorted_vals[-1]
    else:
        if value < sorted_vals[0] or value > sorted_vals[-1]:
            raise Exception("value must (%f) be within the range of values" % value)

    ix = None
    for i in range(len(sorted_vals) - 1):
        # simple snap
        if sorted_vals[i] == value:
            ix = i
        elif sorted_vals[i+1] == value:
            ix = i + 1
        elif sorted_vals[i] < value < sorted_vals[i+1]:
            if direction is 'down':
                ix = i
            elif direction is 'up':
                ix = i + 1
            else:
                dists = [abs(sorted_vals[i] - value), abs(sorted_vals[i+1] - value)]
                ix = i + dists.index(min(dists))
        # satisfy minimum
        if ix is not None:
            if minimum is not None:
                for ii in range(ix,len(sorted_vals)-1):
                    if minimum <= sorted_vals[ii]:
                        ix = ii
                        break
            break
    if ix is None:
        raise Exception('Check if sequence is ordered low to high')
    return sorted_vals[ix]

def find_wet_columns(geo, features, keep_polygons=False):
    """ returns feature_columns, a list of list, each list is col names that the
    corresponding feature goes through.

    If keep_polygons, will return additional list of feature polygons.  Each
    element of the list is a tuple of (feature shape, [column polygons
    intersected]).
    """
    feature_columns = []
    polygons = []
    ### go through features
    for feature in features:
        shp, meta = load_feature(feature['file'])
        cols, pors, polys = get_columns_intersect_polygon(shp, geo,
            threshold=feature['threshold'])
        feature_columns.append(cols)
        # saves (feature polygon, list of polygons of columns intersected)
        polygons.append((shp, polys))
    if keep_polygons:
        return feature_columns, polygons
    else:
        return feature_columns

def line_cross_polygon(line, polygon):
    """ Returns positions of a LineString (eg. a river) where it goes through
    the edges of a Polygon (eg. a mulgraph column).  Note the positions are
    normalised by the length of the LineString.

    If the whole LineString is within the polygon, a tuple of (0.0, 1.0) will be
    returned.

    Requires shapely (LineString and Polygon)
    """
    from shapely.ops import split as sh_split
    p_start = Point(*line.coords[0])
    p_end = Point(*line.coords[-1])
    splitted = sh_split(line, polygon)
    ns = len(splitted)
    if ns >= 3:
        plot_point(Point(*splitted[0].coords[-1]))
        plot_point(Point(*splitted[-1].coords[0]))
        x, y = splitted[-1].coords[0]
        plt.plot(x, y, '+', color='red')
        ndd1 = splitted[0].length / line.length
        ndd2 = 1.0 - (splitted[-1].length / line.length)
        return (ndd1, ndd2)
    elif ns == 2:
        if polygon.contains(p_start):
            ndd = splitted[0].length / line.length
            return (0.0, ndd)
        elif polygon.contains(p_end):
            ndd = 1 - (splitted[-1].length / line.length)
            return (ndd, 1.0)
        else:
            raise Exception
    elif ns == 1:
        # the column does not split the line, ie the line is contained within
        # the column
        if polygon.contains(line):
            return (0.0, 1.0)
        else:
            # line not crossing/inside the column, use project point
            print 'line not crossing/inside the column, use project point'
            nd = line.project(polygon.centroid, normalized=True)
            return (nd, nd)
    else:
        raise Exception

def plot_features(geo, features, polygons=None,
                  column_names=[], column_texts={}):
    """
    polygons is expected to be a list, each element of the list is a tuple of
    (feature shape, [column polygons intersected]). column_names is a list of
    column names to label in plot. column_texts is a dictionary keyed by column
    names, the text will be printed.  Similar to column_names, but can be used
    to print anything.
    """
    geo.layer_plot(plt=plt,
                   column_names=column_names,
                   linecolour='b',
                   title=geo.filename,
                   wells=True, well_names=False)
    for i, feature in enumerate(features):
        if polygons:
            shp, col_polys = polygons[i]
            # plot feature and label
            patch = PolygonPatch(shp, fc=GRAY, ec=GRAY, alpha=0.3, zorder=2)
            plt.gca().add_patch(patch)
            plt.gca().text(shp.centroid.x, shp.centroid.y, feature['name'],
                           verticalalignment='center',
                           horizontalalignment='center',
                           color=GRAY)
            # plot columns that feature intersects
            for p in col_polys:
                patch = PolygonPatch(p, fc=BLUE, ec=BLUE, alpha=0.3, zorder=2)
                plt.gca().add_patch(patch)
    for c in geo.columnlist:
        if c.name in column_texts:
            plt.gca().text(c.centre[0], c.centre[1], column_texts[c.name],
                verticalalignment='center', color=BLACK)

def interp_z(x, z1, z2):
    return z1 + x * (z2-z1)

def plot_polygon(ob, fc='#999999', alpha=0.5):
    patch = PolygonPatch(ob, facecolor=fc, alpha=alpha, zorder=2)
    plt.gca().add_patch(patch)

def plot_point(ob):
    x, y = ob.xy
    plt.plot(x, y, 'o', color='#999999', zorder=1)

def plot_bounds(ob):
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    plt.plot(x, y, 'o', color='#000000', zorder=1)

def plot_line(ob, color='#6699cc', alpha=0.3, linewidth=3):
    x, y = ob.xy
    plt.plot(x, y, color=color, alpha=alpha, linewidth=linewidth,
             solid_capstyle='round', zorder=2)

def plot_col_surface(cname, d1, d2, z, surf, bottom):
    p = Polygon([(d1,bottom),(d1,surf),(d2,surf),(d2,bottom)])
    patch = PolygonPatch(p, facecolor='#999999', alpha=0.5, zorder=2)
    plt.gca().add_patch(patch)

    p = Polygon([(d1,surf),(d1,z),(d2,z),(d2,surf)])
    patch = PolygonPatch(p, edgecolor='blue', alpha=0.2, zorder=1)
    plt.gca().add_patch(patch)

    plt.text((d1+d2)*0.5, z, "%s\n%.2f" % (cname, z),
             horizontalalignment='center',
             verticalalignment='top')

def process_wet_columns(geo, features, feature_columns, overwrite={}):
    """ snaps wet columns for waiwera and calculate the water depths
    """
    column_mins = {}
    wet_cols = {}
    layer_tops = [lay.top for lay in geo.layerlist][::-1]
    for ii,(feature, cols) in enumerate(zip(features, feature_columns)):
        if feature['type'] == 'river':
            for col in cols:
                fitted = geo.column[col].surface
                # new_z = snap(layer_tops, fitted - feature['depth'], 'down')
                # liq_h = fitted - new_z  # wet atm depth
                # nb_min = fitted
                new_z = snap(layer_tops, fitted - feature['depth'], 'down')
                if col in overwrite:
                    new_z = overwrite[col]
                    del overwrite[col]
                    print col, '->', new_z
                liq_h = feature['depth']  # wet atm depth
                nb_min = new_z
                for nb_col in geo.column[col].neighbour:
                    if nb_col.name in column_mins:
                        column_mins[nb_col.name] = max(nb_min, column_mins[nb_col.name])
                    else:
                        column_mins[nb_col.name] = nb_min
                geo.column[col].surface = new_z
                wet_cols[col] = (liq_h, feature['name'])
        elif feature['type'] == 'lake':
            for col in cols:
                fitted = geo.column[col].surface
                new_z = snap(layer_tops, fitted, 'down')
                # special check to ensure lake bottom below elevation
                while new_z >= feature['elevation']:
                    iz = layer_tops.index(new_z)
                    new_z = layer_tops[iz - 1]
                if col in overwrite:
                    new_z = overwrite[col]
                    del overwrite[col]
                    print col, '->', new_z
                liq_h = feature['elevation'] - new_z
                nb_min = feature['elevation']
                for nb_col in geo.column[col].neighbour:
                    if nb_col.name in column_mins:
                        column_mins[nb_col.name] = max(nb_min, column_mins[nb_col.name])
                    else:
                        column_mins[nb_col.name] = nb_min
                geo.column[col].surface = new_z
                wet_cols[col] = (liq_h, feature['name'])
        elif feature['type'] == 'cascade':
            print '+++', feature['name']
            plt.clf()

            cen_line, meta = load_feature(feature['centreline_file'])
            riv_poly, meta = load_feature(feature['file'])
            level = feature['elevation'] # a list expected: [start, end]
            min_depth = feature['min_depth']
            plt.subplot(2,1,1)
            plot_line(cen_line, color='#6699cc', alpha=0.9, linewidth=1)
            plot_polygon(riv_poly, fc='#6699cc', alpha=0.3)
            plt.subplot(2,1,2)
            plt.plot([0.0, cen_line.length], [level[0], level[1]], '-')
            for col in cols:

                fitted = geo.column[col].surface
                col_poly = geo_column_polygon(geo, col)
                plt.subplot(2,1,1)
                nd1, nd2 = line_cross_polygon(cen_line, col_poly)

                plt.subplot(2,1,1)
                plot_polygon(col_poly, fc='#999999', alpha=0.1)
                plt.axis('equal')

                z1 = interp_z(nd1, level[0], level[1])
                z2 = interp_z(nd2, level[0], level[1])
                water_top = (z1 + z2) / 2.0
                new_z = snap(layer_tops, water_top - min_depth, 'down')
                if col in overwrite:
                    new_z = overwrite[col]
                    del overwrite[col]
                    print col, '->', new_z
                liq_h = water_top - new_z  # wet atm depth
                nb_min = water_top
                for nb_col in geo.column[col].neighbour:
                    if nb_col.name in column_mins:
                        column_mins[nb_col.name] = max(nb_min, column_mins[nb_col.name])
                    else:
                        column_mins[nb_col.name] = nb_min
                geo.column[col].surface = new_z
                wet_cols[col] = (liq_h, feature['name'])

                plt.subplot(2,1,2)
                bottom = layer_tops[layer_tops.index(new_z) - 1]
                plot_col_surface(col, cen_line.length*nd1, cen_line.length*nd2,
                                 water_top, new_z, bottom)
                plt.ylim(bottom, None)
                print "Column '%s' from river %.2f -- %.2f, %.2f" % (col, nd1, nd2, water_top),
                print new_z, bottom


            plt.gcf().set_size_inches(10., 10.)
            basename = os.path.splitext(os.path.basename(geo.filename))[0]
            plt.savefig('%s.png' % (basename + '_cascade_wet_features_' + str(ii)))
            # plt.show()

        else:
            raise Exception('process_wet_columns(): Unrecognised feature type %s' % feature['type'])
    for col in geo.columnlist:
        # we columns, adjust, need to take care of wet_cols
        if col.name in wet_cols:
            continue
        # dry columns
        if col.name in overwrite:
            new_z = overwrite[col.name]
            print col.name, '-->', new_z
        else:
            if col.name in column_mins:
                new_z = snap(layer_tops, col.surface, None, column_mins[col.name])
            else:
                new_z = snap(layer_tops, col.surface, None, None, allow_outside=True)
        col.surface = new_z
    return wet_cols, column_mins

if __name__ == '__main__':
    g_base_file = 'gtmp_ay2017_05_5a'

    geo = mulgrid(g_base_file + '.dat')
    with open('wet_surface_features.json', 'r') as f:
        data = json.load(f)
        features = data['features']

    ### find columns that features intersects
    feature_columns, polygons = find_wet_columns(geo, features, True)
    cols_fit_min = list(set([c for cols in feature_columns for c in cols]))

    ### plot
    plot_features(geo, features, polygons, column_names=cols_fit_min)
    plt.gcf().set_size_inches(37., 21.)
    plt.savefig('%s_features.png' % geo.filename)
    plt.show()

    ### fit surface
    dem = np.fromfile('Wairakei_incl_Lake.xyz',sep=" ")
    nrow = np.size(dem) / 3
    dem1= dem.reshape(nrow,3)
    geo.fit_surface(dem1, alpha =0.01, beta =0.01,
        min_columns=cols_fit_min,
        grid_boundary=True)

    geo.write(g_base_file + '_fit.dat')
    with open('_features_columns.json', 'w') as f:
        json.dump(feature_columns, f, indent=4, sort_keys=True)

    ### snapping to layer boundaries (Waiwera needs this)
    wet_cols, column_mins = process_wet_columns(geo, features, feature_columns)

    geo.write(g_base_file + '_snap.dat')
    with open('_wet_column_depth.json', 'w') as f:
        json.dump(wet_cols, f, indent=4)

    ### plotting the results of snapping
    v_top_bycol, v_i_bycol = {}, {}
    for col in geo.columnlist:
        if col.name not in wet_cols.keys() + column_mins.keys():
            continue
        v_i_bycol[col.name] = col.num_layers
        v_top_bycol[col.name] = col.surface
        if col.name in wet_cols:
            v_top_bycol[col.name] += wet_cols[col.name][0]
        # convert to string
        v_i_bycol[col.name] = '%i' % v_i_bycol[col.name]
        v_top_bycol[col.name] = '%.1f' % v_top_bycol[col.name]

    plot_features(geo, features, polygons, column_texts=v_i_bycol)
    plt.gcf().set_size_inches(37., 21.)
    plt.savefig('%s.png' % (geo.filename+'_snap_i'))
    plt.show()

    plot_features(geo, features, polygons, column_texts=v_top_bycol)
    plt.gcf().set_size_inches(37., 21.)
    plt.savefig('%s.png' % (geo.filename+'_snap_surf'))
    plt.show()

