from mulgrids import *
from matplotlib import pyplot as plt
from descartes import PolygonPatch

from shapely.geometry import shape
from shapely.geometry import Polygon

import json
import unittest

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
    if cache_shapely:
        if hasattr(geo, '_column_polygons_'):
            c_polygons = geo._column_polygons_
        else:
            c_polygons = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]
            geo._column_polygons_ = c_polygons
    else:
        c_polygons = [Polygon([n.pos for n in c.node]) for c in geo.columnlist]

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

def process_wet_columns(geo, features, feature_columns, overwrite={}):
    """ snaps wet columns for waiwera and calculate the water depths
    """
    column_mins = {}
    wet_cols = {}
    layer_tops = [lay.top for lay in geo.layerlist][::-1]
    for feature, cols in zip(features, feature_columns):
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

