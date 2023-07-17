import numpy as np


def compute_similariteAire(courbe1, courbe2):
    simi = np.sum(1 / (1 + (courbe1 - courbe2) ** 2))
    simi /= len(courbe1)
    return simi


def compute_Rcorrelation(courbe1, courbe2, verbose=False):
    R = np.corrcoef(courbe2, courbe1)[0, 1]
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
            print("good model")
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
            print("acceptable model")
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


def compute_indicateurComp(courbe1, courbe2, par_R2=0.7, par_FA2=0.8, par_FB=0.3, par_FS=0.05, par_NMSE=0.4,
                           par_simAire=0.95, par_MGmin=0.75, par_MGmax=1.25, v=False):
    mini = np.minimum(np.min(courbe1), np.min(courbe2))
    if mini < 0:
        courbe1 = courbe1 + np.abs(mini)
        courbe2 = courbe2 + np.abs(mini)
    res = {}
    res['R2'] = compute_Rcorrelation(courbe1, courbe2, verbose=v)
    res['FA2'] = compute_erreurFA2(courbe1, courbe2, verbose=v)
    res['FB'] = compute_erreurFB(courbe1, courbe2, verbose=v)
    res['FS'] = compute_erreurFS(courbe1, courbe2, verbose=v)
    res['NMSE'] = compute_erreurNMSE(courbe1, courbe2, verbose=v)
    res['NMSE_O'] = compute_erreurqn(courbe1, courbe2, verbose=v)
    res['MG'] = compute_erreurMG(courbe1, courbe2, verbose=v)
    res['VG'] = compute_erreurVG(courbe1, courbe2, verbose=v)
    res['simAire'] = compute_similariteAire(courbe1, courbe2)
    res['distMax'] = compute_distMaxi(courbe1, courbe2)
    test = True
    test &= res['R2'] >= par_R2
    test &= res['FA2'] >= par_FA2
    test &= np.abs(res['FB']) < par_FB
    test &= np.abs(res['FS']) < par_FS
    test &= res['NMSE'] < par_NMSE
    test &= res['simAire'] > par_simAire
    test &= res['MG'] >= par_MGmin
    test &= res['VG'] >= par_MGmin
    test &= res['VG'] <= par_MGmax
    totT = sum([res[x] >= y for x, y in zip(['R2', 'FA2', 'FB', 'FS', 'NMSE', 'simAire', 'MG', 'VG'],
                                            [par_R2, par_FA2, par_FB, par_FS, par_NMSE, par_simAire, par_MGmin,
                                             par_MGmin])])
    totF = 9 - totT
    if test:
        res['all'] = 1
    else:
        res['all'] = -1
    res['Nombre de conditions vraies'] = totT
    res['Nombre de conditions fausses'] = totF
    return res
