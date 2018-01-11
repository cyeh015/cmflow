from mulgrids import *
from t2grids import *
from t2data import t2generator
from t2data_json import t2data_export_json as t2data

import os.path

def basename_from_geo_filename(geo_filename):
    basename = os.path.splitext(geo_filename)[0]
    if basename.startswith('g'):
        basename = basename[1:]
    return basename

def create_basic_t2data(geo, parameter={}, multi={}, others={}):
    dat = t2data()
    simul = 'AUTOUGH2.2'
    dat.title = ''
    dat.parameter.update(
        {
            'max_timesteps': 999,
            'tstop': 1.e14,
            'const_timestep': 1.e7,
            'print_interval': 999,
            'print_level': 3,
            'gravity': 9.81,
            'default_incons': [101.3e3, 10.0, 0.0],
        })
    # dat.start = True
    # Set MOPs:
    dat.parameter['option'][1] = 1
    dat.parameter['option'][16] = 5
    # set EOS
    dat.multi = {
        'eos': 'EW',
        'num_components': 1,
        'num_equations': 2,
        'num_phases': 2,
        'num_secondary_parameters': 6,
    }
    # add grid, from geo
    dat.grid = t2grid().fromgeo(geo)

    ##### update with user settings #####
    dat.parameter.update(parameter)
    dat.multi.update(multi)
    eos = dat.multi['eos']
    #####################################

    # convert_...() will fail if empty filename
    dat.filename = basename_from_geo_filename(geo.filename) + '.dat'
    dat.convert_to_AUTOUGH2(warn=False, MP=False, simulator=simul, eos=eos)

    return dat


if __name__ == '__main__':
    geo = mulgrid('gtmp_ay2017_05_5a_snap.dat')

    parameter = {
        'default_incons': [101.3e3, 10.0, 0.0],
    }
    multi = {
        'eos': 'EWAV',
        'num_components': 2,
        'num_equations': 3,
        'num_phases': 2,
        'num_secondary_parameters': 6,
    }

    ### create dat from geo
    dat = create_basic_t2data(geo)
    filename = basename_from_geo_filename(geo.filename)
    dat.write(filename + '.dat')

    ### create incon
    dep = dat.parameter['default_incons']
    inc = dat.grid.incons(tuple(dep))
    inc.write(filename + '.incon')

    ### create geners, put into dat
    # additional things like GENERs etc.
    # In my model.json I have things like { "HeatFlux": {"Default": 0.06} }
    model_settings = {
        'HeatFlux': {
            'Default': 0.08
        }
    }
    # add base heat geners
    if 'HeatFlux' in model_settings:
        if 'Default' in model_settings['HeatFlux']:
            heatflow = model_settings['HeatFlux']['Default']
            layername=geo.layerlist[-1].name
            for col in geo.columnlist:
                bname = geo.block_name(layername,col.name)
                gname = geo.block_name('99', col.name)
                gx = col.area * heatflow
                gen=t2generator(name=gname, block=bname, type='HEAT', gx=gx)
                dat.add_generator(gen)

    ### convert to waiwera input (exo and json)
    dat.write_exodus_json(
        geo,
        indent = 2,
        atmos_volume = 1.e25,
        incons = inc,
        eos = None,
        bdy_incons = None,
        mesh_coords = 'xyz')
