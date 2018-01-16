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

def update_rocktype_bycopy(dat, blk_names, to_rocktype, convention='++***'):
    """
        this copies the rocktype onto blocks in blk_names, if need new rocktype
        names, part of the name can be preserved, the rest copied.
        to_rocktype has to be a PyTOUGH 'rocktype' object
    """
    for b in blk_names:
        r_name = ''
        for i,c in enumerate(convention):
            if c == '+':
                r_name = r_name + dat.grid.block[b].rocktype.name[i]
            else:
                r_name = r_name + to_rocktype.name[i]
        if r_name not in dat.grid.rocktype:
            from copy import deepcopy
            new_rock = deepcopy(to_rocktype)
            new_rock.name = r_name

            # THIS IS ONLY FOR EMILY AND ME, DANGEROUS!!!
            for i,c in enumerate(convention[2:5]):
                if c == '+':
                    new_rock.permeability[i] = dat.grid.block[b].rocktype.permeability[i]

            dat.grid.add_rocktype(new_rock)
            print '      new rocktype added: ', new_rock.name
            dat.grid.block[b].rocktype = new_rock
        else:
            dat.grid.block[b].rocktype = dat.grid.rocktype[r_name]

def update_block_geology(dat, blk_name, rock_name):
    """ only updates the block's rocktype name.  rock_name is a 5 chars string,
    can contain '+' character to indicate preservin that part of name.  If the
    final rocktype name does not exist in dat.grid, it will be created by
    copying the current rocktype. """
    def merge_name(orig, new):
        """ both should be 5 chars long, and new can contain '+' """
        final = ''
        for i in range(5):
            if new[i] <> '+':
                final += new[i]
            else:
                final += orig[i]
        return final

    orig_rock_name = dat.grid.block[blk_name].rocktype.name
    new_rock_name = merge_name(orig_rock_name, rock_name)

    if new_rock_name in dat.grid.rocktype:
        dat.grid.block[blk_name].rocktype = dat.grid.rocktype[new_rock_name]
    else:
        from copy import deepcopy
        new_rock = deepcopy(dat.grid.rocktype[orig_rock_name])
        new_rock.name = new_rock_name
        dat.grid.add_rocktype(new_rock)
        dat.grid.block[blk_name].rocktype = new_rock
        print '      new rocktype added: ', new_rock_name

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
