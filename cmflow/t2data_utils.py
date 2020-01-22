from mulgrids import *
from t2grids import *
from t2data import t2generator
from t2data import *
# from t2data_json import t2data_export_json as t2data

import os.path
import json
import re

def basename_from_geo_filename(geo_filename):
    basename = os.path.splitext(geo_filename)[0]
    if basename.startswith('g'):
        basename = basename[1:]
    return basename

def modify_wellname(mod,original):
    """ modifying well name according to rules: any '*' characters in mod
    will be keep original name, otherwise overwritten by mod. """
    newname = []
    for i in xrange(5):
        if mod[i] == '*':
            newname.append(original[i])
        else:
            newname.append(mod[i])
    return "".join(newname)

def gaussian(pos, centre, peak, sigma):
    """ Returns value of 2D Gaussian function:
      f(x,y) = A exp(- ((x-x0)^2+(y-y0)^2)/(2 sigma^2))
    pos is [x, y]; centre is [x0, y0]; peak is A; sigma is standard deviation. """
    from math import exp, pow
    e = (pow((pos[0]-centre[0]),2.0)+pow((pos[1]-centre[1]),2.0))/(2.0 * pow(sigma, 2.0))
    return peak * exp(- e)

def add_heat_geners(geo, dat, fconfig,
                    newwell_label='***99', skip_if_exist='***98'):
    """ adds HEAT geners into dat according to settings in fconfig file """
    # heat flow unit is W/m2
    with open(fconfig, 'r') as f:
        data = json.load(f)

    zone_cols = {}
    for zn,zpoly in data["ZonePolygon"].iteritems():
        if zn not in data["HeatFlux"].keys():
            print "skip unused zone %s" % zn
            continue
        zone_cols[zn] = geo.columns_in_polygon(zpoly)

    layername = geo.layerlist[-1].name
    total_gx = 0.0
    for col in geo.columnlist:
        heatfluxs = []
        # minimum of default value
        for zone in data["HeatFlux"].keys():
            if zone == "Default":
                heatfluxs.append(data["HeatFlux"]["Default"])
            elif zone in zone_cols:
                if col in zone_cols[zone]:
                    heatfluxs.append(data["HeatFlux"][zone])
            elif zone in data["GaussianNodes"]:
                heatfluxs.append(gaussian(
                    col.centroid,
                    data["GaussianNodes"][zone]['centre'],
                    data["HeatFlux"][zone],
                    data["GaussianNodes"][zone]['size']))

        gx = col.area * max(heatfluxs)
        total_gx = total_gx + gx
        blockname = geo.block_name(layername, col.name)
        wellname = modify_wellname(newwell_label, blockname)
        gen=t2generator(name=wellname,block=blockname,type='HEAT',gx=gx)
        # need to check if mass already exist, if does, skip adding
        skip_well = modify_wellname(skip_if_exist, blockname)
        if (blockname, skip_well) not in dat.generator:
            dat.add_generator(gen)
    return dat, total_gx

def add_mass_geners(geo, dat, fconfig,
                    newwell_label='***99', only_blocks=None, only_rocks=None):
    with open(fconfig, 'r') as f:
        data = json.load(f)
    zone_cols = {}
    for zn,zpoly in data["ZonePolygon"].iteritems():
        if zn not in data["Upflow"].keys():
            print "skip unused zone %s" % zn
            continue
        zpoly = [np.array(zz) for zz in zpoly]
        zone_cols[zn] = geo.columns_in_polygon(zpoly)
    layername = geo.layerlist[-1].name
    total_gx = 0.0
    for name,area in data["Upflow"].iteritems():
        blocks = [geo.block_name(layername,col.name) for col in zone_cols[name]]
        if only_blocks is not None:
            pat = re.compile(only_blocks)
            blocks = [b for b in blocks if re.match(pat, b)]
        if only_rocks is not None:
            pat = re.compile(only_rocks)
            blocks = [b for b in blocks if re.match(pat, dat.grid.block[b].rocktype.name)]
        for b in blocks:
            ex = area["Enthalpy"]
            if "Total" in area:
                gx = area["Total"] / float(len(blocks))
            elif "Flux" in area:
                gx = col.area * area["Flux"]
            elif "Amount" in area:
                gx = area["Amount"]
            else:
                raise Exception("Options other than Total, Flux, or Amount is not suppoerted.")
            gn = modify_wellname(newwell_label, b)
            dat.add_generator(t2generator(name=gn, block=b, type='MASS', gx=gx, ex=ex))
            total_gx += gx
    return dat, total_gx

def enthalpy(t_or_p,ph='liq'):
    """ Return enthalpy kJ/kg ('liq', 'vap', or 'dif') of water at specified
        temperature (<=500.0 in degC) or pressu (>500.0 in Pa) """
    import t2thermo
    def enth(t,p,f):
        d,u = f(t,p)
        return u + p/d
    def hlhs((t,p)):
        return enth(t,p,t2thermo.cowat), enth(t,p,t2thermo.supst)
    def sat_tp(t_or_p):
        if t_or_p > 500.0:
            return t2thermo.tsat(t_or_p), t_or_p
        else:
            return t_or_p, t2thermo.sat(t_or_p)
    (hl,hs) = hlhs(sat_tp(t_or_p))
    # xxxx
    return {'liq': hl,'vap': hs,'dif': hs-hl}[ph]

def add_rain_geners(geo, dat, config):
    import t2thermo
    if isinstance(config, str):
        with open(config, 'r') as f:
            config = json.load(f)

    if 'IncludeColumnsOnly' in config and len(config['IncludeColumnsOnly']) != 0:
        cols = config['IncludeColumnsOnly']
    elif 'ExludeColumns' in config and len(config['ExludeColumns']) != 0:
        cols = [c.name for c in geo.columnlist if c.name not in config['ExludeColumns']]
    elif ('ExludeColumns' not in config or len(config['ExludeColumns']) == 0) and ('IncludeColumnsOnly' not in config or len(config['IncludeColumnsOnly']) == 0):
        cols = [c.name for c in geo.columnlist]
    else:
        raise Exception('Error: Both [ExludeColumns] and [IncludeColumnsOnly] used, this is confusing.')

    rain_temp = config['Rain Temperature']
    newlabel = config['NewGenerLabel']
    infiltration = config['Infiltration']
    annualrain = config['AnnualRainFall mm/yr']
    (rain_density,u) = t2thermo.cowat(rain_temp,101325.0)
    rain_enth = enthalpy(rain_temp)
    mmyr2kgs = rain_density/1000.0/365.25/24.0/60.0/60.0

    total_rain = 0.0
    for c in cols:
        col = geo.column[c]
        layername = geo.layerlist[geo.num_layers-col.num_layers].name
        blockname = geo.block_name(layername, col.name)
        genname = modify_wellname(newlabel, blockname)

        rain = col.area * annualrain * infiltration * mmyr2kgs
        total_rain = total_rain + rain

        gen=t2generator(name=genname, block=blockname, type='MASS',
                        gx=rain, ex=rain_enth)
        dat.add_generator(gen)
    return dat, total_rain

def create_rain_geners(geo, config):
    """ Creates two lists of geners, one for NS one for PR.  NS one uses
    constant AnnualRainFall, while PR uses TimedRainFall.  If TimedRainFall does
    not exist in the config, PR will be the same as NS, constant MSS.
    """
    import t2thermo
    if isinstance(config, str):
        with open(config, 'r') as f:
            config = json.load(f)

    if 'IncludeColumnsOnly' in config and len(config['IncludeColumnsOnly']) != 0:
        cols = config['IncludeColumnsOnly']
    elif 'ExludeColumns' in config and len(config['ExludeColumns']) != 0:
        cols = [c.name for c in geo.columnlist if c.name not in config['ExludeColumns']]
    elif ('ExludeColumns' not in config or len(config['ExludeColumns']) == 0) and ('IncludeColumnsOnly' not in config or len(config['IncludeColumnsOnly']) == 0):
        cols = [c.name for c in geo.columnlist]
    else:
        raise Exception('Error: Both [ExludeColumns] and [IncludeColumnsOnly] used, this is confusing.')

    rain_temp = config['Rain Temperature']
    newlabel = config['NewGenerLabel']
    infiltration = config['Infiltration']
    annualrain = config['AnnualRainFall mm/yr']
    (rain_density,u) = t2thermo.cowat(rain_temp,101325.0)
    rain_enth = enthalpy(rain_temp)
    mmyr2kgs = rain_density/1000.0/365.25/24.0/60.0/60.0

    # rain history
    do_pr = False
    if 'TimeOffset yr' in config:
        offset = ['TimeOffset yr']
    else:
        offset = 0.0
    if 'TimedRainFall yr,mm/yr' in config:
        times, rains = zip(*config['TimedRainFall yr,mm/yr'])
        times, rains = np.array(times), np.array(rains)
        times = (times - offset) * 60.*60.*24.*365.25
        enths = np.ones_like(rains)
        enths = enths * rain_enth
        do_pr = True

    total_rain = 0.0
    ns_gs, pr_gs = [], []
    for c in cols:
        col = geo.column[c]
        layername = geo.layerlist[geo.num_layers-col.num_layers].name
        blockname = geo.block_name(layername, col.name)
        genname = modify_wellname(newlabel, blockname)

        # NS
        rain = col.area * annualrain * infiltration * mmyr2kgs
        total_rain = total_rain + rain
        gen=t2generator(name=genname, block=blockname, type='MASS',
                        gx=rain, ex=rain_enth)
        ns_gs.append(gen)

        # PR
        if do_pr:
            rates = col.area * rains * infiltration * mmyr2kgs
            gen=t2generator(name=genname, block=blockname, type='MASS',
                            gx=rain, ex=rain_enth,
                            time=times, rate=rates, enthalpy=list(enths),
                            ltab=len(times), itab='1')
            pr_gs.append(gen)
        else:
            pr_gs.append(gen)

    return ns_gs, pr_gs



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
    return new_rock_name

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
