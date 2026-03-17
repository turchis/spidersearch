import numpy as np
from astropy.timeseries import LombScargle

#These functions give as output the periodograms and False-Alarm-Probability 'FAP1' (see Miles-Paez et al. 2017) using n_bootstrap=1000 for LombScargle, phase dispersion minimization (Stellingwerf 1978) and multi-band phase dispersion minimization (Huijse et al. 2018)

def LombScargle_search(mjd, mag, sigma_mag, P_low, P_up, n_steps, nterms):
    ls = LombScargle(mjd[~np.isnan(mag)], mag[~np.isnan(mag)], sigma_mag[~np.isnan(mag)], nterms=nterms, fit_mean=True)
    freq = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
    freq = 1/freq
    power = ls.power(freq)
    return 1/freq, power

def LombScargle_search_test(mjd, mag_test, P_low, P_up, n_steps, nterms):
    ls = LombScargle(mjd[~np.isnan(mag_test)], mag_test[~np.isnan(mag_test)], nterms=nterms, fit_mean=True)
    freq = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
    freq = 1/freq
    power = ls.power(freq)
    return 1/freq, power

def LombScargle_FAP_astropy(mjd, mag, sigma_mag, P_low, P_up, n_steps, nterms):
    ls = LombScargle(mjd[~np.isnan(mag)], mag[~np.isnan(mag)], sigma_mag[~np.isnan(mag)], nterms=nterms, fit_mean=True)
    freq = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
    freq = 1/freq
    power = ls.power(freq)
    #FAP = ls.false_alarm_probability(power.max(), method='bootstrap', minimum_frequency=freq.min(), maximum_frequency=freq.max())
    FAP = ls.false_alarm_level(0.0000027, method='bootstrap', minimum_frequency=freq.min(), maximum_frequency=freq.max())
    return FAP

def LombScargle_FAP(mjd, mag, sigma_mag, P_low, P_up, n_steps, nterms, n_bootstraps):
    ls = LombScargle(mjd[~np.isnan(mag)], mag[~np.isnan(mag)], sigma_mag[~np.isnan(mag)], nterms=nterms, fit_mean=True)
    freq = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
    freq = 1/freq
    power = ls.power(freq)
    power_max = power.max()

    counter = 0
    for i in range(0,n_bootstraps):
        index = np.random.permutation(np.size(mjd))
        mag_copy = mag[index]
        sigma_mag_copy = sigma_mag[index]
        ls = LombScargle(mjd[~np.isnan(mag_copy)], mag_copy[~np.isnan(mag_copy)], sigma_mag_copy[~np.isnan(mag_copy)], nterms=nterms, fit_mean=True)
        freq = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
        freq = 1/freq
        power = ls.power(freq)
        if power.max()>=power_max:
            counter = counter+1
    FAP = counter/n_bootstraps
    print(counter)
        
    return FAP

def LombScargle_FAP_level(mjd, mag, sigma_mag, P_low, P_up, n_steps, nterms, FAPl_factor,conf_level):
    n_bootstraps = 10000
    n_steps = int(np.round(n_steps/FAPl_factor))
    FAPl = np.zeros(n_steps)
    FAPl_powers = np.zeros((n_bootstraps,n_steps))

    period = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)
    for i in range(0,n_bootstraps):
        index = np.random.permutation(np.size(mjd))
        mag_copy = mag[index]
        sigma_mag_copy = sigma_mag[index]
        ls = LombScargle(mjd[~np.isnan(mag_copy)], mag_copy[~np.isnan(mag_copy)], sigma_mag_copy[~np.isnan(mag_copy)], nterms=nterms, fit_mean=True)
        power = ls.power(1/period)
        FAPl_powers[i,:] = power
    FAPl = np.percentile(FAPl_powers,conf_level,axis=0)
        
    return period, FAPl

def PDM_search(mjd, mag, P_low, P_up, T_0, n_steps, n_bin, min_bin):
    periods = np.empty(n_steps)
    periods[:] = np.nan
    var = np.empty(n_steps)
    var[:] = np.nan
    var_ts = np.var(mag[~np.isnan(mag)])

    P_up = min(P_up,mjd[-1]-mjd[0])

    for i in range(0,n_steps):
        P = P_low + i*(P_up-P_low)/n_steps
        phase = ((mjd-T_0)/P) % 1
        var_tot = 0.
        count_tot = 0
        for j in range(0,n_bin):
            phi = (j+0.5)/n_bin
            m_inrange = []
            count = 0
            for k in range(0,len(phase)):
                if phase[k]>phi-0.5/n_bin and phase[k]<=phi+0.5/n_bin and not np.isnan(mag[k]):
                    m_inrange.append(mag[k])
                    count = count+1
            if count>=min_bin:
                var_tot = var_tot+(count-1)*np.var(m_inrange)
                count_tot = count_tot + count
        periods[i] = P
        var[i] = var_tot/(count_tot-n_bin)/var_ts

    return periods, var

def PDM_FAP(mjd, mag, P_low, P_up, T_0, n_steps, n_bin, min_bin,n_bootstraps):
    periods = np.empty(n_steps)
    periods[:] = np.nan
    var = np.empty(n_steps)
    var[:] = np.nan
    var_ts = np.var(mag[~np.isnan(mag)])

    P_up = min(P_up,mjd[-1]-mjd[0])

    for i in range(0,n_steps):
        P = P_low + i*(P_up-P_low)/n_steps
        phase = ((mjd-T_0)/P) % 1
        var_tot = 0.
        count_tot = 0
        for j in range(0,n_bin):
            phi = (j+0.5)/n_bin
            m_inrange = []
            count = 0
            for k in range(0,len(phase)):
                if phase[k]>phi-0.5/n_bin and phase[k]<=phi+0.5/n_bin and not np.isnan(mag[k]):
                    m_inrange.append(mag[k])
                    count = count+1
            if count>=min_bin:
                var_tot = var_tot+(count-1)*np.var(m_inrange)
                count_tot = count_tot + count
        periods[i] = P
        var[i] = var_tot/(count_tot-n_bin)/var_ts
    var_min = var.min()

    counter = 0
    for i in range(0,n_bootstraps):
        index = np.random.permutation(np.size(mjd))
        mag_copy = mag[index]
        var = np.empty(n_steps)
        var[:] = np.nan
        
        for j in range(0,n_steps):
            P = P_low + j*(P_up-P_low)/n_steps
            phase = ((mjd-T_0)/P) % 1
            var_tot = 0.
            count_tot = 0
            for k in range(0,n_bin):
                phi = (k+0.5)/n_bin
                m_inrange = []
                count = 0
                for l in range(0,len(phase)):
                    if phase[l]>phi-0.5/n_bin and phase[l]<=phi+0.5/n_bin and not np.isnan(mag_copy[l]):
                        m_inrange.append(mag_copy[l])
                        count = count+1
                if count>=min_bin:
                    var_tot = var_tot+(count-1)*np.var(m_inrange)
                    count_tot = count_tot + count
            var[j] = var_tot/(count_tot-n_bin)/var_ts
        if var.min()<=var_min:
            counter = counter+1
    FAP = counter/n_bootstraps
    print(counter)

    return FAP

def PDM_FAP_level(mjd, mag, P_low, P_up, T_0, n_steps, n_bin, min_bin, FAPl_factor, conf_level):
    n_bootstraps = 10000
    n_steps = int(np.round(n_steps/FAPl_factor))
    FAPl = np.zeros(n_steps)
    FAPl_var = np.zeros((n_bootstraps,n_steps))

    var_ts = np.var(mag[~np.isnan(mag)])
    P_up = min(P_up,mjd[-1]-mjd[0])
    period = np.linspace(P_low,min(P_up,mjd[-1]-mjd[0]),n_steps)

    for i in range(0,n_bootstraps):
        index = np.random.permutation(np.size(mjd))
        mag_copy = mag[index]
        var = np.empty(n_steps)
        var[:] = np.nan
        
        for j in range(0,n_steps):
            P = P_low + j*(P_up-P_low)/n_steps
            phase = ((mjd-T_0)/P) % 1
            var_tot = 0.
            count_tot = 0
            for k in range(0,n_bin):
                phi = (k+0.5)/n_bin
                m_inrange = []
                count = 0
                for l in range(0,len(phase)):
                    if phase[l]>phi-0.5/n_bin and phase[l]<=phi+0.5/n_bin and not np.isnan(mag_copy[l]):
                        m_inrange.append(mag_copy[l])
                        count = count+1
                if count>=min_bin:
                    var_tot = var_tot+(count-1)*np.var(m_inrange)
                    count_tot = count_tot + count
            var[j] = var_tot/(count_tot-n_bin)/var_ts
        FAPl_var[i,:] = var
    FAPl = np.percentile(FAPl_var,100-conf_level,axis=0)

    return period, FAPl

def PDM_bisearch(mjd_a, mag_a, mjd_b, mag_b, P_low, P_up, T_0, n_steps, n_bin, min_bin):
    periods = np.empty(n_steps)
    periods[:] = np.nan
    var = np.empty(n_steps)
    var[:] = np.nan
    var_ts_a = np.var(mag_a[~np.isnan(mag_a)])
    var_ts_b = np.var(mag_b[~np.isnan(mag_b)])

    P_up = min(P_up,max(mjd_a[-1]-mjd_a[0],mjd_b[-1]-mjd_b[0]))

    for i in range(0,n_steps):
        P = P_low + i*(P_up-P_low)/n_steps
        phase_a = ((mjd_a-T_0)/P) % 1
        var_tot_a = 0.
        count_tot_a = 0
        phase_b = ((mjd_b-T_0)/P) % 1
        var_tot_b = 0.
        count_tot_b = 0
        for j in range(0,n_bin):
            phi = (j+0.5)/n_bin
            
            m_inrange_a = []
            count_a = 0
            for k in range(0,len(phase_a)):
                if phase_a[k]>phi-0.5/n_bin and phase_a[k]<=phi+0.5/n_bin and not np.isnan(mag_a[k]):
                    m_inrange_a.append(mag_a[k])
                    count_a = count_a+1
            if count_a>=min_bin:
                var_tot_a = var_tot_a+(count_a-1)*np.var(m_inrange_a)
                count_tot_a = count_tot_a + count_a

            m_inrange_b = []
            count_b = 0
            for k in range(0,len(phase_b)):
                if phase_b[k]>phi-0.5/n_bin and phase_b[k]<=phi+0.5/n_bin and not np.isnan(mag_b[k]):
                    m_inrange_b.append(mag_b[k])
                    count_b = count_b+1
            if count_b>=min_bin:
                var_tot_b = var_tot_b+(count_b-1)*np.var(m_inrange_b)
                count_tot_b = count_tot_b + count_b
                
        periods[i] = P
        var[i] = var_tot_a/(count_tot_a-n_bin)/var_ts_a + var_tot_b/(count_tot_b-n_bin)/var_ts_b

    return periods, var

def PDM_trisearch(mjd_a, mag_a, mjd_b, mag_b, mjd_c, mag_c, P_low, P_up, T_0, n_steps, n_bin, min_bin):
    periods = np.empty(n_steps)
    periods[:] = np.nan
    var = np.empty(n_steps)
    var[:] = np.nan
    var_ts_a = np.var(mag_a[~np.isnan(mag_a)])
    var_ts_b = np.var(mag_b[~np.isnan(mag_b)])
    var_ts_c = np.var(mag_c[~np.isnan(mag_c)])

    P_up = min(P_up,max(mjd_a[-1]-mjd_a[0],mjd_b[-1]-mjd_b[0],mjd_c[-1]-mjd_c[0]))

    for i in range(0,n_steps):
        P = P_low + i*(P_up-P_low)/n_steps
        phase_a = ((mjd_a-T_0)/P) % 1
        var_tot_a = 0.
        count_tot_a = 0
        phase_b = ((mjd_b-T_0)/P) % 1
        var_tot_b = 0.
        count_tot_b = 0
        phase_c = ((mjd_c-T_0)/P) % 1
        var_tot_c = 0.
        count_tot_c = 0
        for j in range(0,n_bin):
            phi = (j+0.5)/n_bin
            
            m_inrange_a = []
            count_a = 0
            for k in range(0,len(phase_a)):
                if phase_a[k]>phi-0.5/n_bin and phase_a[k]<=phi+0.5/n_bin and not np.isnan(mag_a[k]):
                    m_inrange_a.append(mag_a[k])
                    count_a = count_a+1
            if count_a>=min_bin:
                var_tot_a = var_tot_a+(count_a-1)*np.var(m_inrange_a)
                count_tot_a = count_tot_a + count_a

            m_inrange_b = []
            count_b = 0
            for k in range(0,len(phase_b)):
                if phase_b[k]>phi-0.5/n_bin and phase_b[k]<=phi+0.5/n_bin and not np.isnan(mag_b[k]):
                    m_inrange_b.append(mag_b[k])
                    count_b = count_b+1
            if count_b>=min_bin:
                var_tot_b = var_tot_b+(count_b-1)*np.var(m_inrange_b)
                count_tot_b = count_tot_b + count_b

            m_inrange_c = []
            count_c = 0
            for k in range(0,len(phase_c)):
                if phase_c[k]>phi-0.5/n_bin and phase_c[k]<=phi+0.5/n_bin and not np.isnan(mag_c[k]):
                    m_inrange_c.append(mag_c[k])
                    count_c = count_c+1
            if count_c>=min_bin:
                var_tot_c = var_tot_c+(count_c-1)*np.var(m_inrange_c)
                count_tot_c = count_tot_c + count_c
                
        periods[i] = P
        var[i] = var_tot_a/(count_tot_a-n_bin)/var_ts_a + var_tot_b/(count_tot_b-n_bin)/var_ts_b + var_tot_c/(count_tot_c-n_bin)/var_ts_c

    return periods, var




