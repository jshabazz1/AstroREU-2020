import mastURL
from astropy.io import fits
from astroquery.mast import Catalogs


def gettimeflux_1800(tic_id, sector):
    url = mastURL.gen_url(tic_id, sector)
    fits.getdata(url, ext=1)
    with fits.open(url, mode="readonly") as hdulist:
        time = hdulist[1].data['TIME']
        flux = hdulist[1].data['FLUX_RAW']
    return time, flux

def gettimeflux_120(tic_id, sector):
    manifest = mastURL.dload_url(tic_id, sector)
    url = manifest[1][0]
    fits.getdata(url, ext=1)
    fits.info(url)
    with fits.open(url, mode="readonly") as hdulist:
        time = hdulist[1].data['TIME']
        flux = hdulist[1].data['SAP_FLUX']
    return time, flux

def get_catalog_data(tic_id):
    catalogTIC = Catalogs.query_criteria(ID=tic_id, catalog='Tic', provenance_name='TASOC',)
    return catalogTIC['Tmag'], catalogTIC['Teff'], catalogTIC['logg'],  

def get_header_data(tic_id, sector):
    url = mastURL.gen_url(tic_id, sector)
    with fits.open(url, mode="readonly") as hdulist:
        CCD = hdulist[0].header['CCD']
        Camera = hdulist[0].header['Camera']
    return CCD, Camera

