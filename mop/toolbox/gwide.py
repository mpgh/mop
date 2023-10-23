"""
TODO: add overlay of 2MASS image:
https://irsa.ipac.caltech.edu/applications/2MASS/IM/docs/siahelp.html
"""

from collections import OrderedDict
import warnings
import numpy as np
from matplotlib import pyplot as plt

from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import vstack

import astropy.utils.exceptions
import astroquery.exceptions

warnings.filterwarnings("ignore")

warnings.simplefilter('ignore', category=warnings.WarningMessage)
warnings.simplefilter('ignore', category=astropy.utils.exceptions.AstropyWarning)
warnings.simplefilter('ignore', category=astroquery.exceptions.InputWarning)
warnings.simplefilter('ignore', category=astroquery.exceptions.InvalidQueryError)

def AOstrehl(mag, dist, seeing=0.8, wfs='visible'):
    """
    mag: guide star magnitude in band (G vor visible, K for IR)
    dist: distance from AO star, in arcsec
    seeing: value of seeing in arcsec (at 500nm)
    wfs: wavefront sensro type, 'visible' (default) or 'IR'
    """
    if wfs.lower()=='visible':
        # http://www.eso.org/sci/facilities/paranal/telescopes/vlti/documents/VLT-MAN-ESO-15000-4552_v111.pdf
        # figure 5
        s = np.interp(mag, [-10,   6,   9,    11, 13.5, 14,  14.5, 15, 15.5, 100], 
                           [0.4, 0.4, 0.37, 0.35, 0.31, 0.2, 0.15, 0.1,   0,   0])

    elif wfs.lower()=='ir':
        # http://www.eso.org/sci/facilities/paranal/telescopes/vlti/documents/VLT-MAN-ESO-15000-4552_v111.pdf
        # figure 7
        s = np.interp(mag, [-10,   7,   8, 9, 10, 10.5, 100], 
                           [0.6, 0.6, 0.55, 0.4, 0.1, 0, 0])
    else:
        assert False, "unknown 'wfs', should be 'visibile' or 'IR'"

    # http://www.eso.org/sci/facilities/paranal/telescopes/vlti/documents/VLT-MAN-ESO-15000-4552_v111.pdf
    # figure 7
    s *= np.interp(dist, [0, 5,   10, 15, 20, 28, 30, 35, 45], 
                         [1, 0.8, 0.6, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05]) 
    return s

def FTvisloss(dist, seeing=0.6, theta0=None, D=820):
    """
    dist: distance from FT star, in arcsec
    seeing: value of seeing in arcsec (at 500nm) default is 0.6
    theta0: isopistonic angle in arcsec (at 500nm) default is None
        will be estimated based on seeing statistics in Paranal
    D: telescope diameter, in cm (default 820)
    """
    if theta0 is None:
        # -- isolplanetic angle measured in Paranal (K band), from Gwide commissioning
        # https://www.aanda.org/articles/aa/pdf/2022/09/aa43941-22.pdf
        # table 3
        #c = np.polyfit(np.log([0.52, 0.62, 0.76]), np.log([3.01, 2.48, 1.96]), 1)
        c = [-1.13111414,  0.36411506]
        theta0 = np.exp(np.polyval(c, np.log(seeing)))
    # -- r0 in cm at 500nm based on seeing
    r0 = 0.98*0.5e-6/(seeing*np.pi/180/3600)*1e2
    #print('theta0=%.1f", r0=%.1fcm at 500nm'%(theta0, r0))

    wl = 2.2 # wavelength in um for Gwide
    r0 *=(wl/0.5)**(6/5)
    theta0 *=(wl/0.5)**(6/5)

    #print('theta0=%.1f", r0=%.1fcm at %.1fum'%(theta0, r0, wl))
    sig_p = 0.12*np.pi**(1/3)*wl*(D/r0)**(-1/6)*dist/theta0
    Vloss = np.exp(-2*np.pi**2/wl**2*sig_p**2)
    return Vloss

def colorPercent(value, thresholds=None, norm=100):
    cols = ['\033[45m', '\033[41m', '\033[43m', '\033[42m', '\033[46m', '\033[44m']
    #cols = ['\033[41m', '\033[43m', '\033[42m', '\033[46m', '\033[44m']
    #cols = ['\033[35m', '\033[31m', '\033[33m', '\033[32m', '\033[36m', '\033[34m']
    
    if thresholds is None:
        thresholds = norm*np.linspace(0,1,len(cols)+1)[1:]
    for i,v in enumerate(thresholds):
        if value<=v:
            return cols[i]+'%2.0f'%min(value,99)+'\033[0m'
    return cols[-1]+'%2.0f'%min(value,99)+'\033[0m'

def findOffFromGSC(SC, gmax=14.5, kmax=10.5, distmax=30, fig=None, 
                   FTAOsame=False, allIn=True, debug=False):
    """
    SC: ('hh:mm:ss', 'dd:mm:ss', 'name') or 'name' (for Simbad)
    gmax: max Gaia mag for AO star 
    kmax: max K mag for fringe tracker
    distmax: max distance, in arcsec
    fig: figure for plot (None for no plot)
    FTAOsame: force AO and FT star to be the same
    allIn: force all stars (SC, FT, AO) to be within distmax from one another 
    """
    if debug:
        print('gmax=%f'%gmax, 'kmax=%f'%kmax, 'distmax=%f'%distmax)
    cunits = (u.hourangle, u.deg)

    # -- see https://cdsarc.cds.unistra.fr/viz-bin/cat/I/353
    catalog = 'I/353/gsc242'

    if not fig is None:
        plt.close(fig)
        plt.figure(fig, figsize=(5.5,5))
        ax = plt.subplot(111, aspect='equal')

    ### Create a SkyCoord object for the target object, querying SIMBAD
    ### for the coordinates if necessary
    #print('\n\033[94m', SC, end=' ')
    if type(SC)==str:
        # -- find in Simbad as source
        name = SC
        sim = Simbad.query_object(SC)
        coor = SkyCoord(ra=sim['RA'].value.data[0], dec=sim['DEC'].value.data[0],
                    frame='icrs', unit=(u.hourangle, u.deg))
    else:
        if (type(SC)==tuple or type(SC)==list) and len(SC)>2:
            name = SC[2]
        else:
            name = str(SC)
        # -- assumes coords (ra, dec)
        try:
            coor = SkyCoord("%f %f"%SC[:2], frame='icrs', unit=cunits)
        except:
            coor = SkyCoord("%s %s"%SC[:2], frame='icrs', unit=cunits)
        if not fig is None:
            plt.title(name)
        name = coor

    print('\033[0m')
    x0 = coor.ra.to(u.deg).value
    y0 = coor.dec.to(u.deg).value
    if not fig is None:
        plt.plot(0, 0, '+k', ms=12, alpha=0.5, label='SC', linewidth=1)
        t = np.linspace(0, 2*np.pi, 100)

    ### All plotting code
    Rs = list(10*np.arange(1, int(distmax/10)))
    if not distmax in Rs:
        Rs.append(distmax)
    if not fig is None:
        for r in Rs:
            plt.plot(r*np.cos(t), r*np.sin(t),
                    '--' if r==distmax else ':',
                    color='r' if r==distmax else '0.3',
                    alpha=0.3)
            plt.text(r*0.71, r*0.71, '%.1f"'%r, rotation=45,
                        fontsize=9 if r==distmax else 7,
                        va='center', ha='center', alpha=0.3)
            plt.text(-r*0.71, -r*0.71, '%.1f"'%r, rotation=45,
                        fontsize=9 if r==distmax else 7,
                        va='center', ha='center', alpha=0.3)
            plt.text(-r*0.71, r*0.71, '%.1f"'%r, rotation=-45,
                        fontsize=9 if r==distmax else 7,
                        va='center', ha='center', alpha=0.3)
            plt.text(r*0.71, -r*0.71, '%.1f"'%r, rotation=-45,
                        fontsize=9 if r==distmax else 7,
                        va='center', ha='center', alpha=0.3)
        plt.xlabel(r'E $\leftarrow$ $\delta$RA (")')
        plt.ylabel(r'$\delta$Dec $\rightarrow$ N (")')
        ax.invert_xaxis()

    ### Build Vizier query with slightly different filters, depending on whether
    ### the star used for AO is the same as the FT.  If they are the same object,
    ### only one catalog of results is returned.
    if FTAOsame:
        # -- FT and AO star are same object
        gsc = Vizier(columns=['*', '+_r', 'pmRA', 'pmDE', 'plx']).query_region(name,
                                  radius=(distmax)*u.arcsec, catalog=catalog,
                                  column_filters={"gaiaGmag":"<=%f"%(gmax),
                                                  "tmassKsmag":"<=%f"%kmax})
        gsc[0]['Use'] = 'AO+FT'
    else:
        gscAO = Vizier(columns=['*', '+_r', 'pmRA', 'pmDE', 'plx']).query_region(name,
                                  radius=(distmax)*u.arcsec,
                                  catalog=catalog,
                                  column_filters={"gaiaGmag":"<=%f"%gmax})
        if debug and len(gscAO)>0:
            print('\033[37m', '-'*3, 'AO stars', len(gscAO), '-'*20)
            gscAO[0]['GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag'].pprint(show_unit=False)
            print('\033[0m')

        gscFT = Vizier(columns=['*', '+_r', 'pmRA', 'pmDE', 'plx']).query_region(name,
                                  radius=(distmax)*u.arcsec,
                                  catalog=catalog,
                                  column_filters={"tmassKsmag":"<=%f"%kmax})
        if debug and len(gscFT)>0:
            print('\033[37m', '-'*3, 'Fringe Tracker stars', len(gscFT), '-'*20)
            gscFT[0]['GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag'].pprint(show_unit=False)
            print('\033[0m')
        # gscFT = Vizier(columns=['*', '+_r', 'pmRA', 'pmDE', 'plx']).query_region(name,
        #                           radius=(distmax)*u.arcsec,
        #                           catalog=catalog,
        #                           column_filters={"Ksmag":"<=%f"%kmax})
        #print(gscFT[0].keys())
        # -- interpolate if Jmag or Hmag or W1mag but no Ksmag
        gscFTp = Vizier(columns=['*', '+_r', 'pmRA', 'pmDE', 'plx']).query_region(name,
                                  radius=(distmax)*u.arcsec, catalog=catalog)
        noK = np.array([type(m)!=np.float64 for m in gscFTp[0]['Ksmag']])

        ### Handle exception if no Ks magnitudes are available for the catalog objects
        ###
        if sum(noK):
            gscFTp = gscFTp[0][noK]
            if debug:
                print('trying to interpolate K mag for:')
                gscFTp['GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag'].pprint(show_unit=False)

            withJ  = np.array([type(m)==np.float64 for m in gscFTp['Jmag']])
            withH  = np.array([type(m)==np.float64 for m in gscFTp['Hmag']])
            withW1 = np.array([type(m)==np.float32 for m in gscFTp['W1mag']])
            
            # -- extrapolate from J & H
            w = withJ * withH * ~withW1
            Ksmag = gscFTp['Hmag'].value + 0.2*(gscFTp['Jmag'].value-gscFTp['Hmag'].value)
            for i in range(len(w)):
                if w[i]:
                    gscFTp['Ksmag'][i] = np.round(Ksmag[i], 1)
                    gscFTp['Ksmag'].mask[i] = False

            # -- interpolate from H & W1
            w = ~withJ * withH * withW1
            Ksmag = 0.5*gscFTp['Hmag'].value + 0.5*gscFTp['W1mag'].value
            for i in range(len(w)):
                if w[i]:
                    gscFTp['Ksmag'][i] = np.round(Ksmag[i], 1)
                    gscFTp['Ksmag'].mask[i] = False
                    
            # -- interpolate from J & W1
            w = withJ * ~withH * withW1
            Ksmag = 0.2*gscFTp['Jmag'].value + 0.8*gscFTp['W1mag'].value
            for i in range(len(w)):
                if w[i]:
                    # print(gscFTp['GSC2'][i], 'K~%.1f'%Ksmag[i], 'from J & W1')
                    # print(gscFTp['Ksmag'][i], gscFTp['Ksmag'].mask[i], end='->')
                    gscFTp['Ksmag'][i] = np.round(Ksmag[i], 1)
                    gscFTp['Ksmag'].mask[i] = False
                    # print(gscFTp['Ksmag'][i], gscFTp['Ksmag'].mask[i])
                    
            # -- interpolate from JHW1
            w = withJ * ~withH * withW1
            Ksmag = 0.2*gscFTp['Jmag'].value + 0.5*gscFTp['Hmag'].value + 0.3*gscFTp['W1mag'].value
            for i in range(len(w)):
                if w[i]:
                    gscFTp['Ksmag'][i] = np.round(Ksmag[i], 1)
                    gscFTp['Ksmag'].mask[i] = False

            gscFTp = gscFTp[gscFTp['Ksmag']<=kmax]
            
            if len(gscFTp)>0:
                print('\033[35minterpolated Kmag for:')
                gscFTp['GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag'].pprint(show_unit=False)
                print('\033[0m')
        else:
            gscFTp = []

        if debug:
            print(len(gscAO), len(gscFT), len(gscFTp))


        ### Choose the final output to present, depending on which catalog was produced
        #print(len(gscFT[0]), N)
        if len(gscFT)==0 and len(gscFTp)>0:
            gscFT = gscFTp
        elif len(gscFTp)>0:
            gscFT = vstack([gscFT[0], gscFTp])
        elif len(gscFT)>0:
            gscFT = gscFT[0]
        elif len(gscAO)==0:
            print('no guide stars')
            return False

        ### Calculate the angular separations of the available Guide Stars from the target,
        ### and store the evaluation
        if len(gscFT)==0:
            gsc = gscAO[0]
            gsc['allIn'] = np.zeros(len(gsc), dtype=bool)
            print('no FT stars!')
        if len(gscAO)==0 and len(gscFT)==0:
            print('no AO/FT stars')
            return False
        elif len(gscAO)==0:
            gsc = gscFT
            print('no AO stars!')
            gsc['allIn'] = np.zeros(len(gsc), dtype=bool)
        else:
            gscAO = gscAO[0]
            # -- check is AO,FT,SC and are all close to one another
            in3ft = np.zeros(len(gscFT), dtype=bool)
            in3ao = np.zeros(len(gscAO), dtype=bool)
            if debug:
                print('computing distance AO/FT')
            # -- for all FT
            for i,g in enumerate(gscFT):
                if not in3ft[i] and not g['GSC2'] in gscAO['GSC2']:
                    # -- for all AO
                    for j,h in enumerate(gscAO):
                        distance = np.sqrt(((g['RA_ICRS']-h['RA_ICRS'])*
                            np.cos((g['DE_ICRS'])*np.pi/180))**2+
                                 (g['DE_ICRS']-h['DE_ICRS'])**2)*3600
                        if debug:
                            print(g['GSC2'], h['GSC2'], distance)
                        if distance<=distmax:
                            in3ao[j] = True
                            in3ft[i] = True
                else:
                    in3ft[i] = True
                    for j,h in enumerate(gscAO):
                        if h['GSC2']==g['GSC2']:
                            in3ao[j] = True
                            
            gscFT['allIn'] = in3ft
            gscAO['allIn'] = in3ao

            # -- merge the two
            gsc = gscAO
            for g in gscFT:
                if not g['GSC2'] in gsc['GSC2']:
                    gsc.add_row(g)
            gsc.sort('_r')
            gsc = gsc

    ### Select stars brighter than the limiting magnitude for Guide Stars in G and Ks bands
    gsc['AO'] = [bool(g['Gmag']<=gmax) for g in gsc]
    
    # gsc['SC strehl'][gsc['AO']] = 100*AOstrehl(dist=gsc['_r'][gsc['AO']], 
    #                                     mag=gsc['Gmag'][gsc['AO']])
    gsc['FT'] = [bool(g['Ksmag']<=kmax) for g in gsc]

    if debug:
        print('GSC=', len(gsc))
        gsc['_r', 'GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag', 'allIn'].pprint(show_unit=False)

    # -- restrict to all objects within distmax:
    if any(gsc['allIn']):
        gsc = gsc
        if allIn:
            gsc = gsc[gsc['allIn']]
    else:
        print('\033[35m', end='')
        print('Gmag<%.1f'%gmax, 'Kmag<%.1f'%kmax)
        gsc['_r', 'GSC2', 'Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag'].pprint(show_unit=False)
        print('\033[0m', end='')
        gsc = []
        if not fig is None:
            plt.close(fig)

    ### Tabular output, colour-coded to identify the useable guide stars.
    if len(gsc):
        gsc['Gmag'] = np.round(gsc['Gmag'], 3)
        #print(gsc.columns)
        gsc['VLTI'] = ['AT' if g['Gmag']<12.5 and g['Ksmag']<9.5 else 'UT' for g in gsc]
        _c = [SkyCoord('%f %f'%(g['RA_ICRS'], g['DE_ICRS']), frame='icrs', unit=(u.deg, u.deg))
            for g in gsc]
        clean = lambda x: x.replace('d',':').replace('h',':').replace('m',':').replace('s','')
        gsc['RAhms'] = [clean(c.ra.to(u.hourangle).to_string()[:11]) for c in _c]
        gsc['RAhms'].info.unit='hh:mm:ss.s'
        gsc['DECdms'] = [clean(c.dec.to_string()[:12]) for c in _c]
        gsc['RAhms'].info.unit='dd:mm:ss.s'

        gsc['GSC2', '_r','Gmag', 'Jmag', 'Hmag', 'Ksmag', 'W1mag', 'RAhms', 'DECdms',
                'plx', 'pmRA', 'pmDE', 'AO', 'FT'].pprint(show_unit=True, max_width=200)

        print()
        # == performance matrix =================================
        AO = np.arange(len(gsc))[gsc['AO']]
        FT = np.arange(len(gsc))[gsc['FT']]

        maxStrehl = 40

        print(' '*11+'AO stars |',' '.join([gsc['GSC2'][j]+'  |' for j in AO]))        
        SCstrehl = [100*AOstrehl(dist=gsc['_r'][j], mag=gsc['Gmag'][j]) for j in AO]
        print(' '*20+'| '+
                ' '.join(['SCstrehl:'+colorPercent(s, norm=maxStrehl)+' |' for s in SCstrehl]))
        print('FT stars'+' '*12+'| '+
                ' '.join(['G=%4.1f @%2.0f" |'%(gsc['Gmag'][j], gsc['_r'][j]) for j in AO]))
        # -- line delimiter
        delim = lambda : print('-'*20+'+'+('-'*13+'+')*len(SCstrehl))
        delim()
        for i in FT:
            Vloss = 100*FTvisloss(dist=gsc['_r'][i])
            print(gsc['GSC2'][i]+' '+'SCvis:'+colorPercent(Vloss), end=' | ')
            Kmags = []
            colk = []
            dist = []
            cold = []
            for j in AO:
                # compute strehl for SC and FT
                distAOFT = np.sqrt( ((_c[j].ra.to(u.deg).value - 
                                      _c[i].ra.to(u.deg).value)*
                                    np.cos(_c[i].dec.to(u.rad).value))**2+
                                    (_c[j].dec.to(u.deg).value - 
                                     _c[i].dec.to(u.deg).value)**2)
                distAOFT *= 3600
                dist.append(distAOFT)
                # -- check if corrected Kmag is bright enough
                if dist[-1]==0:
                    cold.append('\033[44m')
                elif dist[-1]<distmax/4:
                    cold.append('\033[46m')
                elif dist[-1]<distmax/2:
                    cold.append('\033[42m')
                elif dist[-1]<3*distmax/4:
                    cold.append('\033[43m')
                elif dist[-1]<=distmax:
                    cold.append('\033[41m')
                else:
                    cold.append('\033[45m')

                FTstrehl = 100*AOstrehl(dist=distAOFT, mag=gsc['Gmag'][j])
                # -- Kmag corrected from strehl
                Kmags.append(gsc['Ksmag'][i] - 2.5*np.log10(FTstrehl/maxStrehl))
                # -- check if corrected Kmag is bright enough
                if Kmags[-1]<kmax-2:
                    colk.append('\033[44m')
                elif Kmags[-1]<kmax-1:
                    colk.append('\033[46m')
                elif Kmags[-1]<kmax:
                    colk.append('\033[42m')
                elif Kmags[-1]<kmax+1:
                    colk.append('\033[43m')
                else:
                    colk.append('\033[41m')
                print('FTstrehl:'+colorPercent(FTstrehl, norm=maxStrehl)+' |', end=' ')
            print() # carriage return
            print('K=%4.1f     @%2.0f"'%(gsc['Ksmag'][i], gsc['_r'][i]), ' '*3, '|', 
                ' '.join(['%sK~%4.1f\033[0m %s@%2.0f"\033[0m |'%(colk[k],K,cold[k], dist[k]) for k,K in enumerate(Kmags)]))
            #delim()
        if not fig is None:
            ra_d  = np.array([c.ra.to(u.deg).value  for c in _c])
            dec_d = np.array([c.dec.to(u.deg).value for c in _c])
            gmag = np.array([g['Gmag'] for g in gsc])
            kmag = np.array([g['Ksmag'] for g in gsc])
            name = np.array([g['GSC2'] for g in gsc])

            w = [g['AO'] or g['FT'] for g in gsc]
            for x,y,g,k,r in zip((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600,
                                (dec_d[w]-y0)*3600, gmag[w], kmag[w], name[w]):
                plt.text(x, y+0.5, '%s\nG%.1f K%.1f'%(r,g,k),
                        ha='center', va='bottom', fontsize=8)

            w = [g['AO'] and g['FT'] for g in gsc]
            plt.plot((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600, (dec_d[w]-y0)*3600,
                     '*', color='green', label='AO+FT')

            w = [g['AO'] and not g['FT'] for g in gsc]
            plt.plot((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600, (dec_d[w]-y0)*3600,
                     '*', color='cyan', label='AO only')

            w = [not g['AO'] and g['FT'] for g in gsc]
            plt.plot((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600, (dec_d[w]-y0)*3600,
                     '*', color='orange', label='FT only')

            w = [not g['AO'] and not g['FT'] for g in gsc]
            plt.plot((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600, (dec_d[w]-y0)*3600,
                     '.', color='0.5', alpha=0.4, ms=6)
            for x,y,g,k,r in zip((ra_d[w]-x0)*np.cos(y0*np.pi/180)*3600,
                                 (dec_d[w]-y0)*3600, gmag[w], kmag[w],
                                 name[w]):
                plt.text(x, y+0.5, '%s\nG%.1f K%.1f'%(r,g,k),
                         ha='center', va='bottom', fontsize=5, alpha=0.4)
            plt.legend(fontsize=7)
            plt.tight_layout()
        return True
    else:
        print('no off-axis source satisfying all criteria')
    return False

