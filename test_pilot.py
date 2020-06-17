from astropy.table import Table
from flatwrm import flatwrm
from matplotlib import pyplot as plt
from Time_and_flux import gettimeflux_1800, get_catalog_data, get_header_data
from astroquery.mast import Observations, Catalogs
from get_IDS import get_IDS
import numpy as np
import csv
from TessCut import my_animation





hardcopy = False
outfile = "/Users/jshabazz/work/textfiles/flares_teff_3000_to_3025.csv"
outfile2 = "/Users/jshabazz/work/textfiles/manifest_teff_3000_to_3025.csv"
debug = False
noplot = False
fit_events = True
magnitude = False
flarepoints = 3
sigma = 4
period=0.
degree=0
fwhm = 0.


#ids = get_IDS(0,634)

#guenther_flares = Table.read("/Users/jshabazz/Work/textfiles/guenther_flares.csv",format='ascii')
#ids = guenther_flares['tic_id'][0:50].astype(str)
ids = {'167602025','141212720'}

#Opens manifest and sets headers 
with open(outfile,"w", newline='') as ofile1, open(outfile2, "w") as ofile2:
    fieldnames = ['TIC', 'Sector','Plot 1', 'Plot 2', 'Teff',
                 'Tmag', 'log(g)','CCD', 'Camera', 'Dwarfs', 
                 'Giants', 'Closest', 'Movie 1', 'Movie 2', 'Movie 3', 'Movie 4']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    ofile2.write(",".join(fieldnames)+'\n')
#
#    fieldnames1 = ["this_id",\
#                 "t_start",\
#                 "t_end",\
#                 "t_max",\
#                 "flux_max",\
#                 "raw_integral",\
#                 "fit_amp",\
#                 "fit_fwhm",\
#                 "fit_t_start",\
#                 "fit_t_end",\
#                 "fit_t_max",\
#                 "fit_integral"]
#    ofile1.write(",".join(fieldnames1)+'\n')
    for this_id in ids:
        try:
            target_name = this_id
            radius = 0.2
            catalogTIC = Catalogs.query_object(target_name, radius, catalog = "TIC")
            numObj = "Number of TIC objects within %f deg of %s: %u" % (radius, target_name, len(catalogTIC))
            where_dwarfs = np.where(catalogTIC['lumclass'] == 'DWARF')[0]
            where_giants = np.where(catalogTIC['lumclass'] == 'GIANT')[0]
            dwarfs = "Number of objects classified as 'DWARF' within %f deg of %s: %u" % (radius, target_name, len(where_dwarfs))
            giants = "Number of objects classified as 'GIANT' within %f deg of %s: %u" % (radius, target_name, len(where_giants))
            where_closest = np.argmin(catalogTIC['dstArcSec'])
            closest = "Closest TIC ID to %s: TIC %s, seperation of %f arcsec. and a TESS mag. of %f" % (target_name, catalogTIC['ID'][where_closest], catalogTIC['dstArcSec'][where_closest], catalogTIC['Tmag'][where_closest])

            #sectors_search = Observations.query_criteria(target_name=this_id, provenance_project='TASOC')
            sectors_search = Observations.query_criteria(target_name=this_id, obs_collection="HLSP", filters="TESS",
                                                         t_exptime=[1799, 1801])
            #print('Getting sectors')
            #import pdb; pdb.set_trace()
            sector_length = len(sectors_search)
            if  sector_length !=0:
                sector = sectors_search['sequence_number']
                sector_length = sector_length-1
                

            print('Getting times')
            #import pdb; pdb.set_trace()
            time = gettimeflux_1800(this_id, str(sector[0]))[0]
            print(time)
            flux = gettimeflux_1800(this_id, str(sector[0]))[1]
            Tmag = get_catalog_data(this_id)[0]
            Teff = get_catalog_data(this_id)[1]
            logg = get_catalog_data(this_id)[2]
            CCD = get_header_data(this_id, str(sector[0]))[0]
            Camera = get_header_data(this_id, str(sector[0]))[1]

            if period == 0:
                period = flatwrm.FindPeriod(time, flux, maxper=10,debug=debug)
            window_var = 1.5

            if degree == 0:
                degree = flatwrm.SelectDegree(time, flux, period * window_var, debug=debug)
            istart, istop = flatwrm.FindFlares(time, flux, period,
                                   returnbinary=False,
                                   N3=flarepoints,
                                   degree=degree,
                                   detection_sigma=sigma,
                                   debug=debug)

            flatwrm.GenerateOutput(time, flux, istart, istop, period,this_id, \
                           fit_events=fit_events, degree=degree, debug=debug, outputfile=outfile)

            #save out put in readme text file

            period = 0.
            degree = 0
            fwhm = 0.

            maxFlares = 4
            numFlares = len(istart)
            if numFlares <= 4:
              n_flares_to_plot = numFlares
            else:
             n_flares_to_plot = maxFlares

            
        except:
            continue



        if numFlares != 0 :


            print('Making plots')
            """Displays a light curve with highlighted flare candidates"""
            figs, ax = plt.subplots()
            ax.scatter(time, flux)
            ax.scatter(time[istart], flux[istart], c = 'm')
            ax.scatter(time[istop], flux[istop], c = 'm')
            figs.suptitle("TIC " + this_id + ": Flare candidates")
            for t0, t1 in zip(istart, istop):
                ax.scatter(time[t0:t1+1], flux[t0:t1+1], c='m')
            plt.savefig('/Users/jshabazz/Work/lightcurves/' + this_id + '.jpg', bbox_inches='tight')
            #plt.show()


            """Displays a grid of plots focused in on each individual flare in the light curve"""
            fig, axs = plt.subplots(2, 2, figsize = (11,7))
            fig.suptitle(f'TIC {this_id}')
            for i in range(n_flares_to_plot):
                if i > 3:
                     print(f'Number of flares detected: {numFlares}')
                     break
                t2 = istart[i]
                t3 = istop[i]
                axs[i//2, i%2].scatter(time, flux, c='k')
                axs[i//2, i%2].scatter(time[t2:t3+1], flux[t2:t3+1])
                axs[i//2, i%2].set_xlim(time[t2]-0.5, time[t3]+0.5)
                #axs[i//2, i%2].set_ylim(flux[t2], flux[t3]+50.0)
                axs[i//2, i%2].set_title(f'Flare candidate: {i+1}')
                my_animation(this_id, time[t2], time[t3], i)
            plt.savefig('/Users/jshabazz/Work/lightcurves/' + this_id + '_flares.jpg', bbox_inches='tight')
            #plt.show()
            ofile2.write(this_id+','+\
                            this_id + '.jpg'+','+\
                            this_id + '_flares.jpg'+','+\
                            str(sector[0])+','+\
                            str(Teff[0])+','+\
                            str(Tmag[0])+','+\
                            str(logg[0])+','+\
                            str(CCD)+','+\
                            str(Camera)+','+\
                            dwarfs+','+\
                            giants+','+\
                            closest+','+\
                            f'{this_id}_flareevent0.gif'+','+\
                            f'{this_id}_flareevent1.gif'+','+\
                            f'{this_id}_flareevent2.gif'+','+\
                            f'{this_id}_flareevent3.gif'+'\n')
            
        else:
            figs, ax = plt.subplots()
            ax.scatter(time, flux)
            figs.suptitle("TIC " + this_id + ": lightcurve")
            plt.savefig('Users/jshabazz/Work/lightcurves_noflares/'+ this_id +'.jpg', bbox_inches='tight')
            continue















