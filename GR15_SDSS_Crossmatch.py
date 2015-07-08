from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.ascii import read
from scipy.integrate import quad
import socket

# Special slice plotting
from demo_floating_axes import setup_axes3

# For changing plotting parameters 
from matplotlib import rcParams
rcParams.update({'axes.facecolor': '#ffffff'})

# Temporary imports
import ipdb


# /--- Functions ---/ #
def ra_to_deg(ra):
    pieces = ra.split(' ')
    hrs = pieces[0]
    mins = pieces[1]
    secs = pieces[2]
    return float(hrs)*(360/24) + float(mins)*(15./60.) + float(secs)*(15./3600.)

def dec_to_deg(dec):
    pieces = dec.split(' ')
    #check to see if unicode minus sign is present
    if dec[0] == u'\u2212':
        degs = -1*float(pieces[0][1:])
        mins = -1*float(pieces[1])
        secs = -1*float(pieces[2])
        return degs+mins/60.+secs/3600.
    # if not, proceed normally
    else:
        degs = float(pieces[0])
        mins = float(pieces[1])
        secs = float(pieces[2])
        return degs+mins/60.+secs/3600.

def Hubble(z,Omega_M=0.3,Omega_L=0.7,h=0.7):
    H0 = 100*h*3.24077929e-20 #mks
    return H0*np.sqrt(Omega_M*((1+z)**3) + Omega_L) 

def ComovingDistance(z,Omega_M=0.3,Omega_L=0.7,inunitsof='hubbledist'):
    H0 = Hubble(0,Omega_M,Omega_L)
    dH = 3e8/H0
    I = quad(lambda Z: 1/(np.sqrt(Omega_M*((1+Z)**3) + Omega_L)),0,z)
    if inunitsof == 'hubbledist':
        return I[0]
    else:
        return I[0]*dH
    
def startup():
    '''
    Retrieves all data from Groener et al. 2015 study.
    '''
    # Opening Excel File
    if socket.gethostname() == 'Umbriel':
        wb = open_workbook('/home/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')
    else:
        wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')
    # First sheet
    sheet_names = wb.sheet_names()
    sheet1 = wb.sheet_by_name(sheet_names[0])
    # Get headers and data
    clusters = sheet1.col_values(0)[1:]
    redshift = sheet1.col_values(1)[1:]
    methods = sheet1.col_values(2)[1:]
    c200 = sheet1.col_values(3)[1:]
    c200_plus = sheet1.col_values(4)[1:]
    c200_minus = sheet1.col_values(5)[1:]
    m200 = sheet1.col_values(6)[1:]
    m200_plus = sheet1.col_values(7)[1:]
    m200_minus = sheet1.col_values(8)[1:]
    cvir = sheet1.col_values(9)[1:]
    cvir_plus = sheet1.col_values(10)[1:]
    cvir_minus = sheet1.col_values(11)[1:]
    mvir = sheet1.col_values(12)[1:]
    mvir_plus = sheet1.col_values(13)[1:]
    mvir_minus = sheet1.col_values(14)[1:]
    short_refs = sheet1.col_values(15)[1:]
    orig_convention = sheet1.col_values(16)[1:]
    cosmology = sheet1.col_values(17)[1:]
    ra = sheet1.col_values(-2)[1:]
    dec = sheet1.col_values(-1)[1:]
    return clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology,ra,dec

def clusters_within_region(ra_min,ra_max,dec_min,dec_max,plotregion=False):
    # mandate that min values are smaller than max values
    assert ra_min < ra_max, "Minimum R.A. must be smaller than maximum value."
    assert dec_min < dec_max, "Minimum Dec. must be smaller than maximum value."
    # all unique cluster names, RA/Dec values, and redshifts
    unique_clusters = [i for i in set(clusters)]
    unique_ra = [ra_to_deg(ra[clusters.index(i)]) for i in unique_clusters]
    unique_dec = [dec_to_deg(dec[clusters.index(i)]) for i in unique_clusters]
    unique_z = [redshift[clusters.index(i)] for i in unique_clusters]
    # select within RA/Dec region
    cl_list = [unique_clusters[i] for i in range(len(unique_clusters))
               if unique_ra[i] >= ra_min
               and unique_ra[i] <= ra_max
               and unique_dec[i] >= dec_min
               and unique_dec[i] <= dec_max]
    ra_list = [unique_ra[i] for i in range(len(unique_clusters))
               if unique_ra[i] >= ra_min
               and unique_ra[i] <= ra_max
               and unique_dec[i] >= dec_min
               and unique_dec[i] <= dec_max]
    dec_list = [unique_dec[i] for i in range(len(unique_clusters))
               if unique_ra[i] >= ra_min
               and unique_ra[i] <= ra_max
               and unique_dec[i] >= dec_min
               and unique_dec[i] <= dec_max]
    z_list = [unique_z[i] for i in range(len(unique_clusters))
               if unique_ra[i] >= ra_min
               and unique_ra[i] <= ra_max
               and unique_dec[i] >= dec_min
               and unique_dec[i] <= dec_max]
    if plotregion is True:
        # do some plotting here
        for i in range(len(cl_list)):
            plt.scatter(ra_list[i],dec_list[i])
        plt.xlim(0,360)
        plt.ylim(-90,90)
        plt.show()
    return cl_list,ra_list,dec_list,z_list
        
def startup_sdss():
    if socket.gethostname() == 'Umbriel':
        fh = read("/home/groenera/Desktop/Dropbox/Private/Research/DataFiles/ClusterEnvironment/sdss_data_trim.csv")
    else:
        fh = read("/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ClusterEnvironment/sdss_data_trim.csv")
    sdss_z = fh['z']
    sdss_ra = fh['ra']
    sdss_dec = fh['dec']
    return sdss_z,sdss_ra,sdss_dec

def plot_dec_slice(dec_min, dec_max, withclusters=False, withbounds=False):

    # prelinimary stuff
    # all unique cluster names, RA/Dec values, and redshifts
    gr15_unique_clusters = [i for i in set(clusters)]
    gr15_unique_ra = [ra_to_deg(ra[clusters.index(i)]) for i in gr15_unique_clusters]
    gr15_unique_dec = [dec_to_deg(dec[clusters.index(i)]) for i in gr15_unique_clusters]
    gr15_unique_z = [redshift[clusters.index(i)] for i in gr15_unique_clusters]
    
    # select from dec slice
    print("Selecting galaxies from declination slice...")
    sdss_dist_trim = [ComovingDistance(sdss_z[i]) for i in range(len(sdss_z)) if sdss_dec[i] >= dec_min and sdss_dec[i] <= dec_max]
    sdss_ra_trim = [sdss_ra[i] for i in range(len(sdss_ra)) if sdss_dec[i] >= dec_min and sdss_dec[i] <= dec_max]
    dist_max = max(sdss_dist_trim)

    # if plotting clusters over-top
    ra_min = 100
    ra_max = 270
    if withclusters is True:
        gr15_ra_trim = [gr15_unique_ra[i] for i in range(len(gr15_unique_z))
                        if gr15_unique_ra[i] >= ra_min
                        and gr15_unique_ra[i] <= ra_max
                        and gr15_unique_dec[i] >= dec_min
                        and gr15_unique_dec[i] <= dec_max]
        gr15_dist_trim = [ComovingDistance(gr15_unique_z[i]) for i in range(len(gr15_unique_z))
                          if gr15_unique_ra[i] >= ra_min
                          and gr15_unique_ra[i] <= ra_max
                          and gr15_unique_dec[i] >= dec_min
                          and gr15_unique_dec[i] <= dec_max]
        gr15_z_trim = [gr15_unique_z[i] for i in range(len(gr15_unique_z))
                          if gr15_unique_ra[i] >= ra_min
                          and gr15_unique_ra[i] <= ra_max
                          and gr15_unique_dec[i] >= dec_min
                          and gr15_unique_dec[i] <= dec_max]
        print("Found {} clusters within this declination slice..".format(len(gr15_ra_trim)))

    fig = plt.gcf()
    ax3, aux_ax3 = setup_axes3(fig, 111, ra0=100, ra1=270, cz1=dist_max)
    aux_ax3.scatter(sdss_ra_trim, sdss_dist_trim,marker='.',s=1,color='black',zorder=1)
    if withclusters is True:
        aux_ax3.scatter(gr15_ra_trim, gr15_dist_trim,marker='o',color='red',zorder=2)
    if withbounds is True:
        radius = 0.0033 # 10 h^-1 Mpc 
        for i in range(len(gr15_ra_trim)):
            #circle=plt.Circle((gr15_ra_trim[i],gr15_dist_trim[i]),radius,color='b',fill=False)
            #aux_ax3.add_artist(circle)
            theta_tmp = np.linspace(0,2*np.pi,400)
            ra_list_tmp = [gr15_ra_trim[i]+(radius/2)*(1+gr15_z_trim[i])/(ComovingDistance(gr15_z_trim[i]))*np.cos(j) for j in theta_tmp]
            z_list_tmp = [gr15_dist_trim[i]+radius*np.sin(j) for j in theta_tmp]
            aux_ax3.plot(ra_list_tmp,z_list_tmp,color='b')
            #fig.gca().add_artist(circle)
    plt.show()
    return

# /--- Preliminary Stuff ---/ #

# Get GR15 clusters
print("Loading cluster data from Groener & Goldberg (2015)...")
clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology,ra,dec = startup()        

# Get sdss galaxy data
print("Loading SDSS galaxy data...")
sdss_z,sdss_ra,sdss_dec = startup_sdss()




if __name__ == "__main__":
    dec_min = 10
    dec_max = 12
    #cl_list,ra_list,dec_list,z_list = clusters_within_region(ra_min=100,ra_max=270,dec_min=dec_min,dec_max=dec_max,plotregion=True)
    plot_dec_slice(dec_min,dec_max,withclusters=True,withbounds=True)
    
