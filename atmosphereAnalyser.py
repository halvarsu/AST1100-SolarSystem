import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.constants import c, Boltzmann, m_p, m_n

plt.style.use('seaborn-whitegrid')
noise_data = np.load('data/sigma_noise.npy')
spectrum_data = np.load('data/spectrum_seed64_600nm_3000nm.npy')

wavelengths, spectrum  = spectrum_data.T
noise = noise_data[:,1]

gases = {'O2':{'lines': [630, 690, 760], 'weight':32 },
         'H20':{'lines': [720, 820, 940], 'weight':18 },
         'CO2':{'lines': [1400, 1600], 'weight': 44},
         'CH4':{'lines': [1660, 2200], 'weight':16 },
         'CO':{'lines': [2340], 'weight': 28},
         'N2O':{'lines': [2870], 'weight': 30}}

res = [300,30,30]

def get_mask(wavelengths,center, width):
    mask  = np.abs(wavelengths-center) < width
    return mask

def get_mask2(wavelengths,min_wl, max_wl):
    mask1 = wavelengths < max_wl
    mask2 = wavelengths > min_wl
    return np.logical_and(mask1, mask2)
    
def plot_spectrums():
    for gas in sorted(gases.keys()):
        info = gases[gas]
        lines = info['lines']
        print '--------', gas, lines,'---------'
        fig, axes = plt.subplots(len(lines))
        try:
            axes[0] 
        except TypeError:
            axes = [axes]
        for i,line in enumerate(lines):
            print gas, line
            mask = get_mask(wavelengths,center = line, width = 0.1)
            axes[i].plot( wavelengths[mask],  spectrum[mask])
            axes[i].axvline(x=line, ymin = 0, ymax = 1, c='r',
                                linewidth = 2)
            axes[i].axvline(x=line+0.015, ymin = 0, ymax = 1, c='g')
            axes[i].legend(['Spectrum around %dnm' %line, '%dnm'%line, 
                            'Possible doppler shift from CO2 line'])
            axes[i].set_title(gas+", "+str(line))
        plt.show()


doppler_found =  {'CO2':{'center':1400,'del_l':0.015,'width': 0.005}}
gases['CO2']['doppler found'] = 0.015

co2_1  = 1400
observed_shift = 0.015
v_sat = c* observed_shift/co2_1
print 'Velocity of satelite if CO2-line is real:', v_sat, 'm/s'

def gas_temp(gas):
    k = Boltzmann
    m = gases[gas]['weight'] *m_p
    d_l = doppler_found[gas]['center']*1e-9
    sigma = doppler_found[gas]['width']*1e-9/(4*np.log(2))
    return m*(sigma*c/d_l)**2/(8*k*np.log(2))

def get_sigma_range(name,line,(T_min,T_max),res = 30):
    k = Boltzmann
    m = gases[name]['weight'] *(m_p+m_n)/2
    lmbd = line #gases[name]['lines'][linenum]
    T = np.linspace(T_min, T_max, res)
    sigma = lmbd/c*np.sqrt((k*T/(m)))
    return sigma 
    
def get_F_range(res=30):
    F_min_min = 0.7
    F_min_max = 1# 0.7
    F_min = np.linspace(F_min_min, F_min_max, res)
    return F_min

def get_l_center_range(name, line ,width= None,res=300):
    v_sat_max = 10000 #m/s
    max_shift = line*v_sat_max/c #nm
    if not width:
        width = max_shift
    mask = get_mask(wavelengths,line, width)
    wl = wavelengths[mask]
    lmbd = np.linspace(wl[0],wl[-1], res)
    return lmbd


def f_model(l, F_min, sig, l_cen):
    F_max = 1
    return F_max - (F_max - F_min)*np.exp(-(l-l_cen)**2/(2*sig**2))


def xhi2(l,F_min, sig, l_cen):
    return np.sum(wavelengths[i])


def search_for_light(name = "CO2", line = 1400):
    temp_range = [150, 450]
    gas_names = sorted(gases.keys())
    m = [gases[name]['weight'] for name in gas_names]
    
    F_values =  get_F_range(res=res[1])
    l_c_values = get_l_center_range(name, line,  res=res[0])
    sigma_values = get_sigma_range(name, line, temp_range, res = res[2])

    sigma_max = np.max(sigma_values)
    max_wl = np.max(l_c_values) + 4*sigma_max
    min_wl = np.min(l_c_values) - 4*sigma_max
    mask = get_mask2(wavelengths, min_wl, max_wl)
    lambda_values = wavelengths[mask]
    spectrum_values = spectrum[mask]
    noise_values = noise[mask]

    xhi = np.zeros(res)
    data = np.array([l_c_values, F_values, sigma_values])
    
    for k,sigma in enumerate(sigma_values):
        print k
        for j,F_min in enumerate(F_values):
            for i,l_c in enumerate(l_c_values):
                mask = get_mask2( lambda_values, l_c - 4*sigma_max, 
                                  l_c + 4*sigma_max)
                lambs = lambda_values[mask]
                fmodl = f_model(lambs, F_min, sigma, l_c)
                val=(spectrum_values[mask]-fmodl)**2/noise_values[mask]**2
                xhi[i,j,k] = np.sum(val)
    return xhi, data

def plot(l_c, F_min, sigma, name, line, args):
    mask = get_mask2(wavelengths,l_c - 8*sigma, l_c + 8*sigma)
    lambds = wavelengths[mask]
    fmodl = f_model(lambds, F_min, sigma, l_c)
    val = (spectrum[mask] - fmodl)**2/noise[mask]**2
    ax = plt.subplot()
    legend = ['model','obs','$\lambda$_center = %f'%l_c]
    ax.plot(lambds, fmodl)
    ax.plot(lambds, spectrum[mask])
    if args.val:
        ax.plot(lambds, val/np.max(val))
        legend.append('xhi')
        plt.ylim([0, max(spectrum[mask]) + 0.12])
    if args.noise:
        ax.plot(lambds, noise[mask])
        legend.append('noise')
        plt.ylim([0, max(spectrum[mask]) + 0.12])
    ax.axvline(x = l_c, ymin=0, ymax = 1, c='0.5',linestyle='--')
    ax.scatter(lambds, spectrum[mask])
    #ax.scatter(l_c, F_min)
    plt.legend(legend)
    plt.xlabel('$\lambda$(nm)')
    plt.ylabel('Normalized flux $F(\lambda)$')
    plt.title("%s :%d" %(name,line))
    plt.xlim([lambds[0],lambds[-1]])
    if args.save_fig:
        filename = 'figure/%s_%d.png'%(name,line)
        print "Saving file ", filename
        plt.savefig(filename)
    plt.show()


def find_spectral_lines(args):
    import glob
    folder = 'data/atmosphericData/'
    iteration = '_'+str(args.iteration)
    if args.verbose:
        print "--------- Using iteration %d -------------" %args.iteration
    names = []
    dl = []
    sigmas = []
    F_mins = []
    from my_argmin import my_argmin
    if args.one_gas:
        chosen_gases = {args.one_gas:gases[args.one_gas]}
    else:
        chosen_gases = gases
    for name in chosen_gases.keys():
        for line in gases[name]['lines']:
            if args.verbose:
                print name, line
            xhi_file = folder+'xhi_%s_%d%s' %(name, line, iteration)
            data_file = folder+'data_%s_%d%s' %(name, line, iteration)
            ext = '.npy'
            if glob.glob(xhi_file+ext):
                if args.verbose == 1:
                    print 'loading existing files:'
                    print '\t-',xhi_file+ext
                    print '\t-',data_file+ext
                xhi = np.load(xhi_file+ext)
                data = np.load(data_file+ext)
            else:
                if args.verbose == 1:
                    print 'creating and saving ', xhi_file+ext
                xhi, data = search_for_light(name, line)
                np.save(xhi_file, xhi)
                np.save(data_file,data)
            i,j,k = my_argmin(xhi)
            if args.verbose == 1:
                print i, j, k
            l_c = data[0][i]
            F_min = data[1][j]
            sigma = data[2][k]
            names.append(name+"_"+str(line))
            dl.append(data[0][i] - line)
            F_mins.append(data[1][j])
            sigmas.append(data[2][k])
            if not args.dont_plot:
                plot(l_c, F_min, sigma, name, line, args)

    print " gas              dl    |   F_min   |   sigma"
    for i, name in enumerate(names):
        print "%-10s : %10f | %10f | %10f "% (name, dl[i],F_mins[i],sigmas[i])

if __name__ == "__main__":
    import sys
    try:
        if sys.argv[1] == 'plot':
            plot_spectrums()
            sys.exit()
        elif sys.argv[1] == 'comp':
            xhi = search_for_light()
            if raw_input('Save files?(y/N)') == 'y':
                np.save('data/xhi_original', xhi)
        elif sys.argv[1] == 'analyse':
            raise NotImplementedError, "ERROR"
            xhi = xhiAnalyser()
        else:
            xhi = np.load('data/xhi_values.npy')
    except IndexError:
        xhi = np.load('data/xhi_values.npy')
        print "no arguments supplied"
    from my_argmin import *
    i,j,k = my_argmin(xhi)
    print i,j,k
