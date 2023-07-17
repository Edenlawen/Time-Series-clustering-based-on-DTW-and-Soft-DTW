import numpy as np

def compute_similariteAire(courbe1, courbe2):
    simi = np.sum(1 / (1 + (courbe1 - courbe2) ** 2))
    simi /= len(courbe1)
    return simi

def compute_Rcorrelation(courbe1, courbe2, verbose=False):
    R = abs(np.corrcoef(courbe2, courbe1)[0, 1])
    if verbose:
        if R ** 2 > 0.9 and R < 0.05:
            print("acceptable model")
        else:
            print("not acceptable model R^2<0.9")
    return R ** 2

def compute_distMaxi(courbe1, courbe2):
    dMax = np.max(np.abs(courbe2 - courbe1))
    return dMax

def compute_erreurFA2(courbe1, courbe2, verbose=False):
    ratio = courbe2 / courbe1
    fraction = ratio[(ratio >= 0.5) & (ratio <= 2)]
    FA2 = len(fraction) / len(courbe1)
    if verbose:
        if FA2 > 0.8:
            print("good model FA2")
        else:
            print("important number of different point")
    return FA2

def compute_erreurFB(courbe1, courbe2, verbose=False):
    m1 = np.mean(courbe1)
    m2 = np.mean(courbe2)
    FB = 2 * (m1 - m2) / (m1 + m2)
    if verbose:
        if np.abs(FB) < 0.3:
            print("acceptable FB")
        else:
            print("non acceptable FB")
    return FB

def compute_erreurFS(courbe1, courbe2, verbose=False):
    sd1 = np.std(courbe1)
    sd2 = np.std(courbe2)
    var1 = sd1 ** 2
    var2 = sd2 ** 2
    FS = 2 * (var1 - var2) / (var1 + var2)
    if verbose:
        if np.abs(FS) < 0.5:
            print("acceptable model FS")
        else:
            print("non acceptable FS")
    return FS

def compute_erreurMG(courbe1, courbe2, verbose=False):
    N = len(courbe1)
    courbe1[courbe1 == 0] = 0.000001
    courbe2[courbe2 == 0] = 0.000001
    m1 = np.mean(np.log(courbe1))
    m2 = np.mean(np.log(courbe2))
    res = np.exp(m1 - m2)
    if verbose:
        if res <= 1.25 and res >= 0.75:
            print("acceptable MG error")
        else:
            print("non acceptable MG error")
    return res

def compute_erreurNMSE(courbe1, courbe2, verbose=False):
    N = len(courbe1)
    Dif = np.mean((courbe1 - courbe2) ** 2)
    m1 = np.mean(courbe1)
    m2 = np.mean(courbe2)
    res = Dif / (m1 * m2)
    if verbose:
        if res < 0.4:
            print("acceptable NMSE error")
        else:
            print("non acceptable NMSE error")
    return res

def compute_erreurqn(courbe1, courbe2, verbose=False):
    N = len(courbe1)
    aire1 = np.sum((courbe1 - courbe2) ** 2)
    aire2 = np.sum(courbe1 ** 2)
    simi = np.sqrt(aire1 / aire2)
    if verbose:
        if simi < 0.4:
            print("acceptable NMSE_O error")
        else:
            print("non acceptable NMSE_O error")
    return simi

def compute_erreurVG(courbe1, courbe2, verbose=False):
    N = len(courbe1)
    courbe1[courbe1 == 0] = 0.000001
    courbe2[courbe2 == 0] = 0.000001
    m1 = np.log(courbe1)
    m2 = np.log(courbe2)
    res = np.exp(np.mean((m1 - m2) ** 2))
    if verbose:
        if res <= 1.25 and res >= 0.75:
            print("acceptable VG error")
        else:
            print("non acceptable VG error")
    return res

def compute_indicateurComp(courbe1, courbe2, par_R2=0.7, par_FA2=0.8, par_FB=0.3, par_FS=0.05,
                           par_NMSE=0.4, par_simAire=0.95, par_MGmin=0.75, par_MGmax=1.25, v=False):
    mini = min(np.min(courbe1), np.min(courbe2))
    if mini < 0:
        courbe1 += abs(mini)
        courbe2 += abs(mini)
    res = {}
    res["R2"] = compute_Rcorrelation(courbe1, courbe2, verbose=v)
    res["FA2"] = compute_erreurFA2(courbe1, courbe2, verbose=v)
    res["FB"] = compute_erreurFB(courbe1, courbe2, verbose=v)
    res["FS"] = compute_erreurFS(courbe1, courbe2, verbose=v)
    res["NMSE"] = compute_erreurNMSE(courbe1, courbe2, verbose=v)
    res["NMSE_O"] = compute_erreurqn(courbe1, courbe2, verbose=v)
    res["MG"] = compute_erreurMG(courbe1, courbe2, verbose=v)
    res["VG"] = compute_erreurVG(courbe1, courbe2, verbose=v)
    res["simAire"] = compute_similariteAire(courbe1, courbe2)
    res["distMax"] = compute_distMaxi(courbe1, courbe2)
    
    test = True
    test = test and (res["R2"] >= par_R2) and (res["FA2"] >= par_FA2) and (np.abs(res["FB"]) < par_FB) and \
           (np.abs(res["FS"]) < par_FS) and (res["NMSE"] < par_NMSE)
    test = test and (res["simAire"] > par_simAire) and (res["MG"] >= par_MGmin) and (res["VG"] >= par_MGmin) and \
           (res["VG"] <= par_MGmax)
    
    totT = 0
    totF = 0
    if res["R2"] >= par_R2:
        totT += 1
    else:
        totF += 1
    if res["FA2"] >= par_FA2:
        totT += 1
    else:
        totF += 1
    if np.abs(res["FB"]) < par_FB:
        totT += 1
    else:
        totF += 1
    if np.abs(res["FS"]) < par_FS:
        totT += 1
    else:
        totF += 1
    if res["NMSE"] < par_NMSE:
        totT += 1
    else:
        totF += 1
    if res["simAire"] > par_simAire:
        totT += 1
    else:
        totF += 1
    if res["MG"] >= par_MGmin:
        totT += 1
    else:
        totF += 1
    if res["VG"] >= par_MGmin:
        totT += 1
    else:
        totF += 1
    if res["VG"] <= par_MGmax:
        totT += 1
    else:
        totF += 1
    
    if test:
        res["all"] = 1
    else:
        res["all"] = -1
    
    res["Nombre de conditions vraies"] = totT
    res["Nombre de conditions fausses"] = totF
    
    return res
