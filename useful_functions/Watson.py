import numpy as np
import nibabel as nib
from scipy.special import erfi
from scipy.optimize import least_squares
from multiprocessing import Pool
from os import listdir
import time


# %% Class
class WMTI_Watson:
    def __init__(self, files, mask=None, params='invivo', nodes=2,
                 rand=False, lb=None, ub=None, md_ub=None, outprefix=''):
        self.params = params
        self.lb = lb
        self.ub = ub
        self.md_ub = md_ub
        self.initialize_parameters()
        self.rand = rand
        self.outprefix = outprefix
        self.nodes = nodes
        if isinstance(files, tuple) or isinstance(files, list):
            inputlist = files
            self.md, self.ad, self.rd = nib.load(inputlist[0]), nib.load(inputlist[1]), nib.load(inputlist[2])
            self.mk, self.ak, self.rk = nib.load(inputlist[3]), nib.load(inputlist[4]), nib.load(inputlist[5])
        elif isinstance(files, str):
            inpath = files
            self.md = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'md' in f.lower()][0])
            self.ad = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'ad' in f.lower()][0])
            self.rd = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'rd' in f.lower()][0])
            self.mk = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'mk' in f.lower()][0])
            self.ak = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'ak' in f.lower()][0])
            self.rk = nib.load([inpath+f'/{f}' for f in listdir(inpath) if 'rk' in f.lower()][0])
        # store headers as outputs
        self.md_affine = self.md.affine
        self.md_header = self.md.header
        # load ROI mask
        if mask is None:
            self.mask = np.logical_not(np.isnan(self.md.get_fdata()))
        else:
            self.mask = nib.load(mask).get_fdata()
        if np.logical_not(
                self.md.shape == self.ad.shape == self.rd.shape == self.mk.shape == self.ak.shape == self.rk.shape == self.mask.shape):
            raise ValueError('Inputs shapes are not consistent. Volumes must have the same shape')

    def __version__(self):
        return print('v1-03.03.22')

    def initialize_parameters(self):
        # Fit lower bound for model parameters[f, Da, Depar, Deperp, kappa]
        if self.params == 'invivo':
            self.lb = [0, 0, 0, 0, 0]
            self.ub = [1, 4, 3, 3, 128]  # fit upper bound for model parameters, in vivo
            self.x0 = [0, 2, 2, 2, 7]  # initial guess, in vivo ### [f, Da, Depar, Deperp, kappa]
            self.md_ub = 2.5  # upper bound on md to avoid CSF contamination
        elif self.params == 'exvivo':
            self.lb = [0, 0, 0, 0, 0]
            self.ub = [1, 2, 2, 2, 128]  # fit upper bound for model parameters, ex vivo
            self.x0 = [0.9, 1.6, 1, 0.4, 7]  # initial guess, ex vivo
            self.md_ub = 1.8
        else:
            if np.all([len(self.lb) == 5, len(self.ub) == 5, len(self.x0) == 5]):
                return  # Use input
            else:
                raise ValueError('Lower bound(lb), upper bound(ub) and initialization parameter(params) must have '
                                 'lenght 5 and follow this order: (f, Da, Depar, Deperp, kappa) or specify one of:'
                                 ' "invivo" or "exvivo"')

    def fit(self):
        t = time.time()
        self.f, self.Da, self.Depar, self.Deperp, self.c2 = WMTI_Watson_maps(self.md.get_fdata(), self.ad.get_fdata(),
                                                                             self.rd.get_fdata(), self.mk.get_fdata(),
                                                                             self.ak.get_fdata(), self.rk.get_fdata(),
                                                                             self.mask, lb=self.lb, ub=self.ub,
                                                                             params=self.x0, md_ub=self.md_ub,
                                                                             rand=self.rand, nodes=self.nodes)
        return print('Completed in ', np.round_(time.time() - t, 3), 's')

    def maps(self):
        if hasattr(self, 'f'):
            return self.f, self.Da, self.Depar, self.Deperp, self.c2
        else:
            AttributeError('WMTI-watson has not fitted maps. Please run .fit() to fit')

    def save(self, outpath):
        if hasattr(self, 'f'):
            for out in ['f', 'Da', 'Depar', 'Deperp', 'c2']:
                newimg = nib.Nifti1Image(eval('self.' + out), affine=self.md_affine, header=self.md_header)
                nib.save(newimg, outpath + f'/{self.outprefix}' + out + '.nii.gz')
        else:
            AttributeError('WMTI-watson has not fitted maps. Please run .fit() to fit')
        return


# %% Functions
def wmti_watson_f(x, moments):
    # moments
    D0 = moments[0]
    D2 = moments[1]
    W0 = moments[2]
    W2 = moments[3]
    W4 = moments[4]

    f = x[0]
    Da = x[1]
    Depar = x[2]
    Deperp = x[3]
    k = x[4]

    dawsonf = 0.5 * np.exp(-k) * np.sqrt(np.pi) * erfi(np.sqrt(k))
    p2 = 1 / 4 * (3 / (np.sqrt(k) * dawsonf) - 2 - 3 / k)
    p4 = 1 / (32 * k ** 2) * (105 + 12 * k * (5 + k) + (5 * np.sqrt(k) * (2 * k - 21)) / dawsonf)

    F1 = 3 * D0 - f * Da - (1 - f) * (2 * Deperp + Depar)
    F2 = 3 / 2 * D2 - p2 * (f * Da + (1 - f) * (Depar - Deperp))
    F3 = D2 ** 2 + 5 * D0 ** 2 * (1 + W0 / 3) - f * Da ** 2 - (1 - f) * (
            5 * Deperp ** 2 + (Depar - Deperp) ** 2 + 10 / 3 * Deperp * (Depar - Deperp))
    F4 = 1 / 2 * D2 * (D2 + 7 * D0) + 7 / 12 * W2 * D0 ** 2 - p2 * (
            f * Da ** 2 + (1 - f) * ((Depar - Deperp) ** 2 + 7 / 3 * Deperp * (Depar - Deperp)))
    F5 = 9 * D2 ** 2 / 4 + 35 / 24 * W4 * D0 ** 2 - p4 * (f * Da ** 2 + (1 - f) * (Depar - Deperp) ** 2)

    return np.array([F1, F2, F3, F4, F5])


def normal_fit_wmti_watson(roi, D0, D2, W0, W2, W4, x0, lb, ub):
    # Empty storage
    fx0, fx1, fx2, fx3, fx4 = [], [], [], [], []
    for i in range(roi[roi].flatten().shape[0]):
        print(str(i) + '/' + str(roi[roi].flatten().shape[0] - 1))

        moments = np.array([D0[roi][i], D2[roi][i], W0[roi][i], W2[roi][i], W4[roi][i]])
        F = least_squares(wmti_watson_f, x0=np.array(x0), bounds=(lb, ub), args=[moments],
                          ftol=1e-6)  # 6 on matlab, ftol

        fx0 += [F.x[0]]
        fx1 += [F.x[1]]
        fx2 += [F.x[2]]
        fx3 += [F.x[3]]
        fx4 += [F.x[4]]
        return fx0, fx1, fx2, fx3, fx4


def parfit_wmti_watson(D0, D2, W0, W2, W4, x0, lb, ub):
    moments = [D0, D2, W0, W2, W4]
    F = least_squares(wmti_watson_f, x0=np.array(x0), bounds=(lb, ub), args=[moments],
                      ftol=1e-12)  # 6 on matlab, ftol
    fx0 = [F.x[0]]
    fx1 = [F.x[1]]
    fx2 = [F.x[2]]
    fx3 = [F.x[3]]
    fx4 = [F.x[4]]

    return fx0, fx1, fx2, fx3, fx4


def rand_x0(leng):
    return np.array([np.random.uniform(0.1, 0.9, leng),
                     np.random.uniform(1.5, 2.5, leng),
                     np.random.uniform(0.5, 1.5, leng),
                     np.random.uniform(0.0, 0.5, leng),
                     np.random.uniform(4.0, 16, leng)]).T     #1/3,1


def WMTI_Watson_maps(md, ad, rd, mk, ak, rk, mask=None, lb=[0, 0, 0, 0, 0], ub=[1, 4, 3, 3, 128], md_ub=2.5,
                     params=[0.9, 2.2, 1.6, 0.7, 7], rand=False, nodes=2):
    '''
    # given md, ad, rd, mk, ak, rk(mean, axial, radial diffusivity, mean, axial, radial kurtosis) maps,
    # calculate WM model parameter maps:
    # f(axonal water fraction), Da(axonal diffusivity), Depar, Deperp(extra - axonal)
    # parallel and perpendicular diffusivities), c2(mean cos ^ 2 of the axon)
    # orientation dispersion: c2 = 1 / 3 fully isotropic, c2 = 1 perfectly parallel)
    # c2 is directly related to the Watson distribution concentration parameter
    # kappa(same as in NODDI)

    # All diffusivities in um2 / ms
    # md, ad, rd should also be in um2 / ms, otherwise converted here

    # mask: brain or ROI mask
    # invivo_flag: boolean, flag for in vivo(true) or ex vivo (false)
    * rand: random initialization

    # I.Jelescu, T. Pavan, Nov. 2021
    '''
    print("Rk : ", rk)
    print("Mk : ", mk)
    print("Ak : ", ak)
    print("Rd : ", rd)
    print("Md : ", md)
    print("Ad : ", ad)
 
    # Avoid /0
    small = 0.0000000001
    x0 = params

    #check if rk <0 and replace with 0
    rk[rk < 0] = small
    mk[mk < 0] = small
    ak[ak < 0] = small

    # Check md in um2 / ms and not mm2 / s, if not, convert
    if np.nanmedian(md) < 1e-2:
        md = md * 1e3
        ad = ad * 1e3
        rd = rd * 1e3

    # Make mask #redundant
    if mask is None:
        mask = np.logical_not(np.isnan(md))
 
    # Filter out voxels with unrealistic tensor values
    filt = (md < md_ub) & (rk > 0) & (rk < 10) & (mk > 0) & (mk < 10)
    roi = np.logical_and(mask, filt)  # exclude voxels with unphysical values from calculation

    # Calculate signal moments
    Wpar = ak * (ad / md) ** 2
    Wperp = rk * (rd / md) ** 2

    D0 = md  # mean diffusivity
    D2 = 2 / 3 * (ad - rd)
    W0 = mk
    W2 = 1 / 7 * (3 * Wpar + 5 * mk - 8 * Wperp)
    W4 = 4 / 7 * (Wpar - 3 * mk + 2 * Wperp)

    # Random initialization?
    if rand:
        prep_x0 = rand_x0(D0[roi].shape[0])
    else:
        prep_x0 = np.tile(x0, (D0[roi].shape[0], 1))

    # initialize model parameter maps
    f = np.zeros_like(md)
    Da = f.copy()
    Depar = f.copy()
    Deperp = f.copy()
    kappa = f.copy()

    # Fit
    # Prepare parameters for parallelization
    zipped = list(zip(D0[roi], D2[roi], W0[roi], W2[roi], W4[roi],
                      prep_x0,
                      np.tile(lb, (D0[roi].shape[0], 1)),
                      np.tile(ub, (D0[roi].shape[0], 1))))
    p = Pool(nodes)
    mapped = np.array(p.starmap(parfit_wmti_watson, zipped)).squeeze()
  
    # Store data in place
    f[roi] = mapped[:, 0]
    Da[roi] = mapped[:, 1]
    Depar[roi] = mapped[:, 2]
    Deperp[roi] = mapped[:, 3]
    kappa[roi] = mapped[:, 4]

    # from kappa, calculate c2, the mean cos ^ 2 of the angle between axons and
    # main bundle orientation(an easier metric of orientations dispersion, c2
    # varies between 1 / 3(isotropic) and 1(perfectly parallel axons)
    Fs = np.sqrt(np.pi) / 2 * np.exp(-kappa) * erfi(np.sqrt(kappa))
    c2 = 1 / (2 * np.sqrt(kappa) * Fs + small) - 1 / (2 * kappa + small)

    return f, Da, Depar, Deperp, np.array(c2)