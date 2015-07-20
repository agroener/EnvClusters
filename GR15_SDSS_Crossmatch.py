from xlrd import open_workbook
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from astropy.io.ascii import read
from scipy.integrate import quad
from scipy.stats import pearsonr
import socket
import time

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

def startup_processed():
    if socket.gethostname() == 'Umbriel':
        fh = read("/home/groenera/Desktop/GitHubRepos/EnvClusters/GR15_normalized.csv")
    else:
        fh = read("/Users/groenera/Desktop/GithubRepositories/envclusters/EnvClusters/GR15_normalized.csv")
    mvir = fh['mvir']
    mvir_p = fh['mvir_p']
    mvir_m = fh['mvir_m']
    cvir = fh['cvir']
    cvir_p = fh['cvir_p']
    cvir_m = fh['cvir_m']
    methods = fh['methods']
    z = fh['z']
    cl = fh['cl']
    refs = fh['refs']
    return mvir,mvir_p,mvir_m,cvir,cvir_p,cvir_m,methods,z,cl,refs

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

def plot_dec_slice(dec_min, dec_max, withclusters=False, withbounds=False, justgalsinside=False, skipplot=False):

    # prelinimary stuff
    # all unique cluster names, RA/Dec values, and redshifts
    gr15_unique_clusters = [i for i in set(clusters)]
    gr15_unique_ra = [ra_to_deg(ra[clusters.index(i)]) for i in gr15_unique_clusters]
    gr15_unique_dec = [dec_to_deg(dec[clusters.index(i)]) for i in gr15_unique_clusters]
    gr15_unique_z = [redshift[clusters.index(i)] for i in gr15_unique_clusters]
    
    # select from dec slice
    print("  --> Selecting galaxies from declination slice...")
    sdss_dist_trim = [ComovingDistance(sdss_z[i]) for i in range(len(sdss_z)) if sdss_dec[i] >= dec_min and sdss_dec[i] <= dec_max]
    sdss_ra_trim = [sdss_ra[i] for i in range(len(sdss_ra)) if sdss_dec[i] >= dec_min and sdss_dec[i] <= dec_max]
    dist_max = max(sdss_dist_trim)

    # if plotting clusters over-top
    ra_min = 100
    ra_max = 270
    if withclusters is True:
        gr15_cl_trim = [gr15_unique_clusters[i] for i in range(len(gr15_unique_z))
                        if gr15_unique_ra[i] >= ra_min
                        and gr15_unique_ra[i] <= ra_max
                        and gr15_unique_dec[i] >= dec_min
                        and gr15_unique_dec[i] <= dec_max]
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
        print("  --> Found {} clusters within this declination slice..".format(len(gr15_ra_trim)))

    fig = plt.gcf()
    ax3, aux_ax3 = setup_axes3(fig, 111, ra0=100, ra1=270, cz1=dist_max)
    if justgalsinside is False:
        aux_ax3.scatter(sdss_ra_trim, sdss_dist_trim,marker='.',s=1,color='black',zorder=1)
    if withclusters is True:
        aux_ax3.scatter(gr15_ra_trim, gr15_dist_trim,marker='o',color='red',zorder=2)
    if withbounds is True:
        radius = 0.0033 # 10 h^-1 Mpc (in units of c/H_0)
        theta_0 = 85*(2*np.pi/360)
        circle_list = []
        for i in range(len(gr15_ra_trim)):
            tmp_ra_rads = gr15_ra_trim[i] * (2*np.pi)/360.
            circle=plt.Circle((-1*gr15_dist_trim[i]*np.cos(tmp_ra_rads+theta_0),-1*gr15_dist_trim[i]*np.sin(tmp_ra_rads+theta_0)),radius,color='b',transform=aux_ax3.transData._b,fill=False)
            circle_list.append(circle)
            aux_ax3.add_artist(circle)
    if justgalsinside is True:
        if skipplot is False:
            print("Plotting only sdss galaxies which are within cluster bounds...")
        master_gal = [[] for i in range(len(circle_list))]
        master_cl = [[] for i in range(len(circle_list))]
        for i in range(len(circle_list)):
            master_cl[i].append([gr15_ra_trim[i], gr15_dist_trim[i], gr15_cl_trim[i]])
            for j in range(len(sdss_ra_trim)):
                tmp_ra_rads = gr15_ra_trim[i] * (2*np.pi)/360.
                #ipdb.set_trace()
                #if circle_list[i].contains_point((-1*gr15_dist_trim[i]*np.cos(tmp_ra_rads+theta_0),-1*gr15_dist_trim[i]*np.sin(tmp_ra_rads+theta_0))):
                if sdss_dist_trim[j] <= gr15_dist_trim[i] + radius and sdss_dist_trim[j] >= gr15_dist_trim[i] - radius:
                    if sdss_ra_trim[j] <= gr15_ra_trim[i] + (radius*(1+gr15_z_trim[i])/ComovingDistance(gr15_z_trim[i]))*(360/(2*np.pi)) and sdss_ra_trim[j] >= gr15_ra_trim[i] - (radius*(1+gr15_z_trim[i])/ComovingDistance(gr15_z_trim[i]))*(360/(2*np.pi)):
                        master_gal[i].append([sdss_ra_trim[j], sdss_dist_trim[j]])
                        aux_ax3.scatter(sdss_ra_trim[j], sdss_dist_trim[j],marker='.',s=1,color='black',zorder=1)
    if skipplot is False:
        plt.show()
    return master_gal,master_cl

def measure_angles(master_gal,master_cl):
    n_cl_in_slice = len(master_cl)
    alpha_list = [[] for i in range(len(master_cl))]
    for i in range(n_cl_in_slice):
        theta_c = master_cl[i][0][0]*(2*np.pi/360)
        r_c = master_cl[i][0][1]
        vec_c = [r_c*np.cos(theta_c),r_c*np.sin(theta_c)]
        vec_c_norm = [k/np.linalg.norm(vec_c) for k in vec_c]
        n_gals_per_cl = len(master_gal[i])
        for j in range(n_gals_per_cl):
            tmp_r_g = master_gal[i][j][1]
            tmp_theta_g = master_gal[i][j][0]*(2*np.pi/360)
            tmp_vec_n = [(tmp_r_g*np.cos(tmp_theta_g) - r_c*np.cos(theta_c)),(tmp_r_g*np.sin(tmp_theta_g) - r_c*np.sin(theta_c))]
            tmp_alpha = math.acos(np.dot(tmp_vec_n,vec_c_norm)/np.linalg.norm(tmp_vec_n))
            if tmp_alpha >= np.pi/2:
                ## range centered around zero (-np.pi/2 to np.pi/2)
                #tmp_alpha = tmp_alpha - np.pi
                ## range from zero to np.pi/2
                tmp_alpha = np.pi - tmp_alpha
            alpha_list[i].append(tmp_alpha)
    
    return alpha_list

def mad(data, axis=None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)

def process_alpha_list(alpha_list, master_gal, master_cl, thresh=100):
    gal_nums_per_cl = [len(master_gal[i]) for i in range(len(master_gal))]
    indices_thresh = [i for i in range(len(alpha_list)) if gal_nums_per_cl[i] >= thresh]
    alpha_list_thresh = [alpha_list[i] for i in indices_thresh]
    cl_thresh = [master_cl[i][0][-1] for i in indices_thresh]
    ## Need to rethink using averages/standard deviations for measuring distribution of angles
    # aves and stds
    aves_thresh = [np.average(alpha_list_thresh[i]) for i in range(len(alpha_list_thresh))]
    stds_thresh = [np.std(alpha_list_thresh[i]) for i in range(len(alpha_list_thresh))]
    # meds and mads
    med_thresh = [np.median(alpha_list_thresh[i]) for i in range(len(alpha_list_thresh))]
    mad_thresh = [mad(alpha_list_thresh[i]) for i in range(len(alpha_list_thresh))]
    # percentiles of gals with alpha <= np.pi/4
    perc_thresh= [np.float(len([alpha_list_thresh[j][i] for i in range(len(alpha_list_thresh[j])) if alpha_list_thresh[j][i] <= np.cos(np.pi/4)]))/len(alpha_list_thresh[j]) for j in range(len(alpha_list_thresh))]
    return cl_thresh,aves_thresh,stds_thresh,med_thresh,mad_thresh,perc_thresh 

def get_concs_and_masses(cl_list, method=None):
    if method is None:
        pro_cl_list = pro_cl.data.tolist()
        #Get all measurements
        indices = [[i for i in range(len(pro_cl_list)) if cl_list[j] == pro_cl_list[i]] for j in range(len(cl_list))]
        p_mvir_list = pro_mvir.data.tolist()
        p_mvir = [[p_mvir_list[i] for i in indices[j]] for j in range(len(indices))]
        p_mvir_p_list = pro_mvir_p.data.tolist()
        p_mvir_p = [[p_mvir_p_list[i] for i in indices[j]] for j in range(len(indices))]
        p_cvir_list = pro_cvir.data.tolist()
        p_cvir = [[p_cvir_list[i] for i in indices[j]] for j in range(len(indices))]
        p_cvir_p_list = pro_cvir_p.data.tolist()
        p_cvir_p = [[p_cvir_p_list[i] for i in indices[j]] for j in range(len(indices))]
        p_z_list = pro_z.data.tolist()
        p_z = [[p_z_list[i] for i in indices[j]] for j in range(len(indices))]
        #Coadd measurements
        coadd_p_mvir = []
        coadd_p_mvir_p = []
        coadd_p_y = []
        coadd_p_y_p = []
        for i in range(len(indices)):
            if len(indices[i]) > 1:
                #print("More than one cluster here: {}".format(i))
                # Coadding masses
                mweights = [(1.0/p_mvir_p[i][j]**2) for j in range(len(p_mvir_p[i]))]
                mnew = sum([p_mvir[i][j]*mweights[j] for j in range(len(p_mvir[i]))])/sum(mweights)
                mnew_p = 1.0/np.sqrt(sum(mweights))
                coadd_p_mvir.append(mnew)
                coadd_p_mvir_p.append(mnew_p)
                # Coadding concs*(1+z)
                y = [p_cvir[i][j]*(1+p_z[i][0]) for j in range(len(p_cvir[i]))]
                yerr = [p_cvir_p[i][j]*(1+p_z[i][0]) for j in range(len(p_cvir_p[i]))]
                yweights = [1.0/yerr[j]**2 for j in range(len(yerr))]
                ynew = sum([y[j]*yweights[j] for j in range(len(yweights))])/sum(yweights)
                ynew_p = 1.0/np.sqrt(sum(yweights))
                coadd_p_y.append(ynew)
                coadd_p_y_p.append(ynew_p)
            else:
                #print("Only one cluster here: {}".format(i))
                coadd_p_mvir.append(p_mvir[i][0])
                coadd_p_mvir_p.append(p_mvir_p[i][0])
                y = [p_cvir[i][j]*(1+p_z[i][0]) for j in range(len(p_cvir[i]))]
                yerr = [p_cvir_p[i][j]*(1+p_z[i][0]) for j in range(len(p_cvir_p[i]))]
                coadd_p_y.append(y[0])
                coadd_p_y_p.append(yerr[0])
                        
        return coadd_p_mvir,coadd_p_mvir_p,coadd_p_y,coadd_p_y_p
    else:
        ipdb.set_trace()
        #Select measurements based upon method
        #Coadd measurements

def correlate_all_slices(dec_min=0,dec_max=66):
    dec_step = 2 # in degrees
    dec_list = range(dec_min,dec_max,dec_step)
    out_mvir,out_mvir_p = ([],[])
    out_y,out_y_p = ([],[])
    out_aves_thresh,out_stds_thresh,out_med_thresh,out_mad_thresh,out_perc_thresh = ([],[],[],[],[])
    for dec in dec_list:
        print("  Declination slice: {} to {}".format(dec,dec+dec_step))
        # Procedure: (1) Select SDSS galaxies in slice; (2) measure angles; (3) correlate with cluster mass/conc measurements
        tmp_dec_min = dec
        tmp_dec_max = dec+dec_step
        # (1)
        master_gal,master_cl = plot_dec_slice(tmp_dec_min,tmp_dec_max,withclusters=True,withbounds=True,justgalsinside=True,skipplot=True)
        # (2)
        alpha_list = measure_angles(master_gal,master_cl)
        cl_thresh,aves_thresh,stds_thresh,med_thresh,mad_thresh,perc_thresh=process_alpha_list(alpha_list,master_gal,master_cl)
        # (3)
        coadd_p_mvir,coadd_p_mvir_p,coadd_p_y,coadd_p_y_p=get_concs_and_masses(cl_thresh)
        out_mvir.append(coadd_p_mvir)
        out_mvir_p.append(coadd_p_mvir_p)
        out_y.append(coadd_p_y)
        out_y_p.append(coadd_p_y_p)
        out_aves_thresh.append(aves_thresh)
        out_stds_thresh.append(stds_thresh)
        out_med_thresh.append(med_thresh)
        out_mad_thresh.append(mad_thresh)
        out_perc_thresh.append(perc_thresh)
    # Flatten data into a single list
    out_mvir = [item for sublist in out_mvir for item in sublist]
    out_mvir_p = [item for sublist in out_mvir_p for item in sublist]
    out_y = [item for sublist in out_y for item in sublist]
    out_y_p = [item for sublist in out_y_p for item in sublist]
    out_aves_thresh = [item for sublist in out_aves_thresh for item in sublist]
    out_stds_thresh = [item for sublist in out_stds_thresh for item in sublist]
    out_med_thresh = [item for sublist in out_med_thresh for item in sublist]
    out_mad_thresh = [item for sublist in out_mad_thresh for item in sublist]
    out_perc_thresh = [item for sublist in out_perc_thresh for item in sublist]
    return out_mvir,out_mvir_p,out_y,out_y_p,out_aves_thresh,out_stds_thresh,out_med_thresh,out_mad_thresh,out_perc_thresh

# /--- Preliminary Stuff ---/ #

# Get all GR15 cluster data
print("Loading cluster data from Groener & Goldberg (2015)...")
clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology,ra,dec = startup()        

# Get normalized concs/masses
print("Loading normalized cluster data (post-processed)...")
pro_mvir,pro_mvir_p,pro_mvir_m,pro_cvir,pro_cvir_p,pro_cvir_m,pro_methods,pro_z,pro_cl,pro_refs = startup_processed()

# Get sdss galaxy data
print("Loading SDSS galaxy data...\n")
sdss_z,sdss_ra,sdss_dec = startup_sdss()




if __name__ == "__main__":
    '''
    # Procedure: (1) Select SDSS galaxies in slice; (2) measure angles; (3) correlate with cluster mass/conc measurements
    dec_min = 10
    dec_max = 12
    # (1)
    master_gal,master_cl = plot_dec_slice(dec_min,dec_max,withclusters=True,withbounds=True,justgalsinside=True)
    # (2)
    alpha_list = measure_angles(master_gal,master_cl)
    cl_thresh,aves_thresh,stds_thresh,med_thresh,mad_thresh,perc_thresh=process_alpha_list(alpha_list)
    # (3)
    coadd_p_mvir,coadd_p_mvir_p,coadd_p_y,coadd_p_y_p=get_concs_and_masses(cl_thresh)
    ipdb.set_trace()
    '''
    t0 = time.time()
    out_mvir,out_mvir_p,out_y,out_y_p,out_aves_thresh,out_stds_thresh,out_med_thresh,out_mad_thresh,out_perc_thresh = correlate_all_slices()
    t1 = time.time()
    print("Total Time: {} seconds.".format(str(t1-t0)))

    
    # Do plotting and stats here
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(out_perc_thresh,np.log10(out_mvir),yerr=np.log10(out_mvir_p),color='red',fmt='o')
    axarr[0].set_ylabel(r"$\mathrm{\log_{10} \, M_{vir}/10^{14} M_{\odot}}$",fontsize=18)
    axarr[0].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[0].set_xlim(0,1)
    axarr[0].text(0.1,1.7,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_perc_thresh,out_mvir)[0],3),round(pearsonr(out_perc_thresh,out_mvir)[1],3)))
    axarr[1].errorbar(out_perc_thresh,np.log10(out_y),yerr=np.log10(out_y_p),color='blue',fmt='o')
    axarr[1].set_ylabel(r"$\mathrm{log_{10}\, c_{vir} \, (1+z)}$",fontsize=18)
    axarr[1].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[1].set_xlim(0,1)
    axarr[1].set_xlabel(r"$\mathrm{\frac{N(\alpha < \pi/4)}{N_{tot}}}$",fontsize=18)
    axarr[1].text(0.1,1.5,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_perc_thresh,np.log10(out_y))[0],3),round(pearsonr(out_perc_thresh,np.log10(out_y))[1],3)))
    plt.show()
   
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(out_aves_thresh,np.log10(out_mvir),yerr=np.log10(out_mvir_p),color='red',fmt='o')
    axarr[0].set_ylabel(r"$\mathrm{\log_{10} \, M_{vir}/10^{14} M_{\odot}}$",fontsize=18)
    axarr[0].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[0].set_xlim(0,np.pi/2)
    axarr[0].text(0.1,1.7,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_aves_thresh,out_mvir)[0],3),round(pearsonr(out_aves_thresh,out_mvir)[1],3)))
    axarr[1].errorbar(out_aves_thresh,np.log10(out_y),yerr=np.log10(out_y_p),color='blue',fmt='o')
    axarr[1].set_ylabel(r"$\mathrm{log_{10}\, c_{vir} \, (1+z)}$",fontsize=18)
    axarr[1].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[1].set_xlim(0,np.pi/2)
    axarr[1].set_xlabel(r"$\mathrm{\langle \alpha \rangle}$",fontsize=18)
    axarr[1].text(0.1,1.5,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_aves_thresh,np.log10(out_y))[0],3),round(pearsonr(out_aves_thresh,np.log10(out_y))[1],3)))
    plt.show()

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(out_med_thresh,np.log10(out_mvir),yerr=np.log10(out_mvir_p),color='red',fmt='o')
    axarr[0].set_ylabel(r"$\mathrm{\log_{10} \, M_{vir}/10^{14} M_{\odot}}$",fontsize=18)
    axarr[0].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[0].set_xlim(0,np.pi/2)
    axarr[0].text(1.2,1.7,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_med_thresh,out_mvir)[0],3),round(pearsonr(out_med_thresh,out_mvir)[1],3)))
    axarr[1].errorbar(out_med_thresh,np.log10(out_y),yerr=np.log10(out_y_p),color='blue',fmt='o')
    axarr[1].set_ylabel(r"$\mathrm{log_{10}\, c_{vir} \, (1+z)}$",fontsize=18)
    axarr[1].set_ylim(-0.1,2) #be careful about chopping data points
    axarr[1].set_xlim(0,np.pi/2)
    axarr[1].set_xlabel(r"$\mathrm{median(\alpha)}$",fontsize=18)
    axarr[1].text(1.2,1.7,"pearsonr: {}\np-value:   {}".format(round(pearsonr(out_med_thresh,np.log10(out_y))[0],3),round(pearsonr(out_med_thresh,np.log10(out_y))[1],3)))
    plt.show()
    
    
    ipdb.set_trace()
