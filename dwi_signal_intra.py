import os, argparse, array, tqdm, parse
import numpy as np 
import multiprocessing as mp 

import dipy.reconst.dki as dki
from dipy.core.gradients import gradient_table


# CONSTANTS 

GYRO        = 267.51525e3 # rad / ms / T
PI          = np.pi
PI_2        = 2*PI
S_TO_MS     = 1e3 
M_TO_MUM    = 1e6
M_TO_MM     = 1e3
MM_TO_MUM   = 1e3

K_COMB_IDX = {1: (0, 0, 0, 0),
 2: (1, 0, 0, 0),
 3: (2, 0, 0, 0),
 4: (1, 1, 0, 0),
 6: (2, 1, 0, 0),
 9: (2, 2, 0, 0),
 8: (1, 1, 1, 0),
 12: (2, 1, 1, 0),
 18: (2, 2, 1, 0),
 27: (2, 2, 2, 0),
 16: (1, 1, 1, 1),
 24: (2, 1, 1, 1),
 36: (2, 2, 1, 1),
 54: (2, 2, 2, 1),
 81: (2, 2, 2, 2)}

K_IND_ELE = {1: 0, 16: 1, 81: 2, 2: 3, 3: 4, 8: 5, 24: 6, 27: 7, 54: 8, 4: 9, 9: 10, 36: 11, 6: 12, 12: 13, 18: 14}

W2_COMB_IDX  = {1: (0, 0), 
                2: (0, 1),
                3: (0, 2), 
                4: (1, 1), 
                6: (1, 2), 
                9: (2, 2)}
W2_IND_ELE  = { 1: 0,
                2: 1,
                3: 2, 
                4: 3, 
                6: 4,
                9: 5}

# 1. PROTOCOL 

def give_b(exp_proto, round=True):

    gyro    = GYRO*S_TO_MS
    b       = (exp_proto[:, 3]*exp_proto[:, 5]*gyro)**2 * (exp_proto[:, 4]-exp_proto[:, 5]/3)

    b       = np.around(b, 6)*(1/M_TO_MM**2) if round else b*(1/M_TO_MM**2)    

    return b # s/mm2

def give_bmat_steam(proto):

    delta_d         = proto.delta.reshape(-1, 1, 1)
    tau_m           = proto.tau_m.reshape(-1, 1, 1)
    delta_c         = proto.delta_c.reshape(-1, 1, 1)
    delta_r         = proto.delta_r.reshape(-1, 1, 1)
    bvecs_steam     = proto.gdir
    G               = proto.G
    cG              = proto.cG
    rG              = proto.rG
    dG              = bvecs_steam * G.reshape(-1, 1)
    tau1            = proto.tau1.reshape(-1, 1, 1)
    tau2            = proto.tau2.reshape(-1, 1, 1)

    tdd = tau1 + tau2 + tau_m + 2*delta_c + 2*delta_d/3 + 2*delta_r
    tcc = tau_m + 2*delta_c/3 + 2*delta_r
    trr = tau_m + 2*delta_r/3
    tdc = tau_m + delta_c + 2*delta_r
    tdr = tau_m + delta_r
    tcr = tau_m + delta_r

    gdgd = np.matmul( np.expand_dims(dG, 2), np.expand_dims(dG, 1))
    gcgc = np.matmul( np.expand_dims(cG, 2), np.expand_dims(cG, 1))
    grgr = np.matmul( np.expand_dims(rG, 2), np.expand_dims(rG, 1))

    gdgc = np.matmul( np.expand_dims(dG, 2), np.expand_dims(cG, 1))
    gcgd = np.matmul( np.expand_dims(cG, 2), np.expand_dims(dG, 1))
    
    gdgr = np.matmul( np.expand_dims(dG, 2), np.expand_dims(rG, 1))
    grgd = np.matmul( np.expand_dims(rG, 2), np.expand_dims(dG, 1))
    
    grgc = np.matmul( np.expand_dims(rG, 2), np.expand_dims(cG, 1))
    gcgr = np.matmul( np.expand_dims(cG, 2), np.expand_dims(rG, 1))

    beff = delta_d**2 * tdd * gdgd + delta_c**2 * tcc * gcgc + delta_r**2 * trr * grgr + delta_d * delta_c * tdc * (gdgc + gcgd) + delta_d*delta_r * tdr * (gdgr + grgd) + delta_c*delta_r * tcr * (gcgr + grgc)
    beff *= (GYRO*S_TO_MS/M_TO_MM)**2 # s/mm2    
    return beff 

def give_bmat_steam_txt(proto):

    delta_d         = proto[:, 5].reshape(-1, 1, 1)
    tau_m           = proto[:, 17].reshape(-1, 1, 1)
    delta_c         = proto[:, 7].reshape(-1, 1, 1)
    delta_r         = proto[:, 11].reshape(-1, 1, 1)
    bvecs_steam     = proto[:, :3]
    G               = proto[:, 3].reshape(-1, 1, 1)
    cG              = proto[:, 8:11]
    rG              = proto[:, 12:15]
    dG              = bvecs_steam * G.reshape(-1, 1)
    tau1            = proto[:, 15].reshape((-1, 1, 1))
    tau2            = proto[:, 16].reshape((-1, 1, 1))

    tdd = tau1 + tau2 + tau_m + 2*delta_c + 2*delta_d/3 + 2*delta_r
    tcc = tau_m + 2*delta_c/3 + 2*delta_r
    trr = tau_m + 2*delta_r/3
    tdc = tau_m + delta_c + 2*delta_r
    tdr = tau_m + delta_r
    tcr = tau_m + delta_r

    gdgd = np.matmul( np.expand_dims(dG, 2), np.expand_dims(dG, 1))
    gcgc = np.matmul( np.expand_dims(cG, 2), np.expand_dims(cG, 1))
    grgr = np.matmul( np.expand_dims(rG, 2), np.expand_dims(rG, 1))

    gdgc = np.matmul( np.expand_dims(dG, 2), np.expand_dims(cG, 1))
    gcgd = np.matmul( np.expand_dims(cG, 2), np.expand_dims(dG, 1))
    
    gdgr = np.matmul( np.expand_dims(dG, 2), np.expand_dims(rG, 1))
    grgd = np.matmul( np.expand_dims(rG, 2), np.expand_dims(dG, 1))
    
    grgc = np.matmul( np.expand_dims(rG, 2), np.expand_dims(cG, 1))
    gcgr = np.matmul( np.expand_dims(cG, 2), np.expand_dims(rG, 1))

    beff = delta_d**2 * tdd * gdgd + delta_c**2 * tcc * gcgc + delta_r**2 * trr * grgr + delta_d * delta_c * tdc * (gdgc + gcgd) + delta_d*delta_r * tdr * (gdgr + grgd) + delta_c*delta_r * tcr * (gcgr + grgc)
    beff *= (GYRO*S_TO_MS/M_TO_MM)**2 # s/mm2    
    return beff 

def give_beff_steam(proto):
    bmat = give_bmat_steam(proto)
    beff = np.expand_dims(np.trace(bmat, axis1=1, axis2=2), 0)
    return beff 

def give_beff_steam_txt(proto):
    bmat = give_bmat_steam_txt(proto)
    beff = np.expand_dims(np.trace(bmat, axis1=1, axis2=2), 0)
    return beff 


    
class Protocol():

    def __init__(self, path_proto, orig, seqtype=None):
    
        self.type   = seqtype

        if orig=="text":
            self.from_text(path_proto)
            self.set_b0_mask()

        elif orig=="matlab":
            self.from_matlab(path_proto)
        else:
            pass

        self.b_ord_idx   = np.argsort(self.b.flatten())
        self.set_val_idx() 
        
        return 

    def from_text(self, path_proto):
        proto, seqtype = read_scheme(path_proto)

        if self.type==None:self.type=seqtype

        if self.type=="PGSE":
            # gx gy gz G Delta delta TE
            self.TE     = proto[:, 6].reshape((1, -1)) * S_TO_MS
            self.G      = proto[:, 3].reshape((1, -1))/M_TO_MUM
            self.gdir   = proto[:, :3]
            self.tdiff  = proto[:, 4].reshape((1, -1))* S_TO_MS
            self.delta  = proto[:, 5].reshape((1, -1))* S_TO_MS
            self.b      = give_b(proto, round=False).reshape((1, -1))* (S_TO_MS/MM_TO_MUM**2)
            self.b0_idx = None

        elif self.type=="STEAM":
            # gx gy gz G Delta delta TE delta_c cGx cGy cGz delta_r rGx rGy rGz tau1 tau2 taum

            self.TE     = proto[:, 6].reshape((1, -1)) * S_TO_MS
            self.G      = proto[:, 3].reshape((1, -1))/M_TO_MUM
            self.gdir   = proto[:, :3]
            self.tdiff  = proto[:, 4].reshape((1, -1))* S_TO_MS
            self.delta  = proto[:, 5].reshape((1, -1))* S_TO_MS
            self.b      = give_beff_steam_txt(proto).reshape((1, -1))* (S_TO_MS/MM_TO_MUM**2)

            self.tau_m      = proto[:, 17].reshape((1, -1))* S_TO_MS
            self.delta_c    = proto[:, 7].reshape((1, -1))* S_TO_MS
            self.delta_r    = proto[:, 11].reshape((1, -1))* S_TO_MS
            self.cG         = proto[:, 8:11]/M_TO_MUM
            self.rG         = proto[:, 12:15]/M_TO_MUM
            self.tau1       = proto[:, 15].reshape((1, -1))* S_TO_MS
            self.tau2       = proto[:, 16].reshape((1, -1))* S_TO_MS
            self.b0_idx = None


        return 


    def from_matlab(self, path_proto):

        from scipy.io import loadmat
        
        d               = loadmat(path_proto)["protocol"]
        proto_mat       = d[0][0]

        self.delta      = proto_mat[1].reshape((1, -1))* S_TO_MS
        self.tdiff      = proto_mat[2].reshape((1, -1))* S_TO_MS
        self.tau_m      = proto_mat[3].reshape((1, -1))* S_TO_MS
        self.delta_c    = proto_mat[4].reshape((1, -1))* S_TO_MS
        self.delta_r    = proto_mat[5].reshape((1, -1))* S_TO_MS
        self.cG         = proto_mat[6]/M_TO_MUM
        self.rG         = proto_mat[7]/M_TO_MUM
        self.gdir       = proto_mat[8]
        self.G          = proto_mat[9].reshape((1, -1))/M_TO_MUM
        self.b          = proto_mat[10].reshape((1, -1)) * (S_TO_MS/M_TO_MUM**2)
        self.TE         = proto_mat[12].reshape((1, -1))* S_TO_MS
        self.tau1       = proto_mat[15].reshape((1, -1))* S_TO_MS
        self.tau2       = proto_mat[16].reshape((1, -1))* S_TO_MS
        self.b0_idx     = proto_mat[17].reshape((1, -1)) - 1

        self.set_b0_mask(d["dataset"][0][0]-1)
        return
    
    def set_val_idx(self, S=None):

        self.val_idx = np.full(self.b.size, True)

        if np.shape(S):
            self.val_idx    = (S > 0 ) & (S<1.2)
        return
    
    def set_b0_mask(self, b0_mask=None):

        if np.shape(b0_mask):
            self.b0_mask = np.array(b0_mask)
        else:
            self.b0_mask = None
        return 
    
    def set_b0_idx(self, b0_idx):
        self.b0_idx = np.array(b0_idx).reshape((1, -1))
        return
    
    def normalize_sig(self, S):

        S_normed = S.copy()

        if np.shape(self.b0_mask):
            for b0i in np.unique(self.b0_mask.flatten()):

                idx_    = self.b0_mask.flatten() == b0i
                idx_b0  = np.intersect1d(np.where(idx_)[0], self.b0_idx)

                if idx_.sum():
                    curr_b0 = S[idx_b0].mean()
                    S_normed[idx_] = S[idx_]/curr_b0

        return S_normed

    def reduce(self, idx_):
        self.TE     = self.TE[idx_].reshape((1, -1))
        self.G      = self.G[idx_].reshape((1, -1))
        self.gdir   = self.gdir[idx_.flatten()]
        self.tdiff  = self.tdiff[idx_].reshape((1, -1))
        self.delta  = self.delta[idx_].reshape((1, -1))
        self.b      = self.b[idx_].reshape((1, -1))

        idx_L           = np.arange(idx_.size)
        old_to_new_idx  = idx_L[idx_.flatten()]

        if self.b0_idx:
            new_bo = []
            for o in self.b0_idx.flatten(): 
                ni = np.where(old_to_new_idx == o)[0]
                if np.size(ni):
                    new_bo.append(ni)
            self.b0_idx  = np.array(new_bo).reshape(1, -1)

        try:
            self.tau_m      = self.tau_m[idx_].reshape((1, -1))
            self.delta_c    = self.delta_c[idx_].reshape((1, -1))
            self.delta_r    = self.delta_r[idx_].reshape((1, -1))
            self.cG         = self.cG[idx_.flatten()]
            self.rG         = self.rG[idx_.flatten()]
            self.tau1       = self.tau1[idx_].reshape((1, -1))
            self.tau2       = self.tau2[idx_].reshape((1, -1))
        except:
            pass

        if self.b0_mask:
            self.b0_mask    = self.b0_mask[idx_] if self.b0_mask is not None else None

        return        

    # Save STEAM Proto in .txt format for Steam from MATLAB
    def proto_steam_to_txt(self, path_save):

        """ 
        Conversion STEAM to PGSE
        pad_i               = proto_steam_A.tau1 + proto_steam_A.delta + proto_steam_A.delta_c + proto_steam_A.delta_r
        proto_pgse_A_c.TE   = proto_steam_A.delta + proto_steam_A.tdiff + proto_steam_A.TE  -pad_i*2 
        """ 

        with open(path_save,"w") as f:
            f.write("VERSION: STEAM\n")

            nb_acq = self.G.size
            for i in range(nb_acq):

                l = "{} "*18 + "\n"
                l = l.format(self.gdir[i, 0], self.gdir[i, 1], self.gdir[i, 2], self.G[0,i]*1e6, self.tdiff[0,i]*1e-3, self.delta[0,i]*1e-3, self.TE[0,i]*1e-3, self.delta_c[0, i]*1e-3, self.cG[i, 0]*1e6, self.cG[i, 1]*1e6, self.cG[i, 2]*1e6, self.delta_r[0, i]*1e-3, self.rG[i, 0]*1e6, self.rG[i, 1]*1e6, self.rG[i, 2]*1e6, self.tau1[0,i]*1e-3, self.tau2[0,i]*1e-3, self.tau_m[0,i]*1e-3)
                f.write(l)

    def proto_pgse_to_txt(self, path_save):
        with open(path_save,"w") as f:
            f.write("VERSION: PGSE\n")

            nb_acq = self.G.size
            for i in range(nb_acq):
                l = "{} "*7 + "\n"
                l = l.format(self.gdir[i, 0], self.gdir[i, 1], self.gdir[i, 2], self.G[0,i]*1e6, self.tdiff[0,i]*1e-3, self.delta[0,i]*1e-3, self.TE[0,i]*1e-3)
                f.write(l)



# 2. READ FILES

def read_scheme(path_scheme):

    with open(path_scheme) as f: lines = f.readlines()

    exp_type = lines[0].split('VERSION: ')[-1].split('\n')[0].split(' ')[0]
    if (exp_type == 'STEJSKALTANNER') |  (exp_type == '1') |  ('PGSE' in exp_type ):
        # gx gy gz G Delta delta TE
        exp_type    = "PGSE"
        exp_proto   = np.zeros((len(lines[1:]), 7))

        for i,line in enumerate(lines[1:]):            
            exp_proto[i] = [float(a) for a in line.split()]

    elif "STEAM" in exp_type:
        # gx gy gz G Delta delta TE delta_c cGx cGy cGz delta_r rGx rGy rGz tau1 tau2 taum

        exp_type    = "STEAM"
        exp_proto   = np.zeros((len(lines[1:]), 18))

        for i,line in enumerate(lines[1:]):            
            exp_proto[i] = [float(a) for a in line.split()]

    return exp_proto, exp_type

def read_info(path_file):
    vox_str     = '( {:g} {:g} {:g} )\n'
    gen_num_str = '{:g}\n'
    gen_bool_str= ' {}\n'
    diff_str    = '{:g}e'
    t_str       = '{:g} ms\n'


    conf = {}
    conf["permeability"] = 0
    with open(path_file) as f: lines = f.readlines()

    for idx_,l in enumerate(lines):
        if 'Voxel'in l:
            conf["vox_min"] = np.array([a for a in parse.parse(vox_str, lines[idx_+1])])*1e-3
            conf["vox_max"] = np.array([a for a in parse.parse(vox_str, lines[idx_+2])])*1e-3

        if 'Number of particles:' in l:
            conf["N"] = int(parse.parse(gen_num_str, l.split('-')[-1])[0])
        
        if 'Number of steps' in l:
            conf["T"] = int(parse.parse(gen_num_str, l.split('-')[-1])[0])
        
        if 'Number of cores' in l:
            conf["num_process"] = int(parse.parse(gen_num_str, l.split('-')[-1])[0])
        
        if 'Diffusivity' in l:
            conf["diffusivity"] = parse.parse(diff_str, l.split('-')[-2])[0] * 1e-9
        
        if 'Particle dynamics duration' in l:
            conf["duration"] = parse.parse(t_str, l.split('-')[-1])[0] * 1e-3

        if 'Permeability' in l:
            conf["permeability"] = 0 if l.split('-')[-1] == " false" else 0
        
        if 'Write to binary' in l:
            tmp = parse.parse(gen_bool_str, l.split('-')[-1])[0]
            conf["bin"] = True if tmp=='true' else False

        if 'Write to txt' in l:
            tmp = parse.parse(gen_bool_str, l.split('-')[-1])[0]
            conf["txt"] = True if tmp=='true' else False

        if 'Walkers initial position' in l:
            tmp = parse.parse(gen_bool_str, l.split('-')[-1])[0]
            conf["ini_walkers_pos"] = tmp

    if not("ini_walkers_pos" in conf.keys()):conf["ini_walkers_pos"] = "delta" 

    return conf

def read_spheres(path_exp, scale_i = 1e-3):
    
    lines = []

    with open(path_exp, 'r') as f:
        
        scale = 1
        for l in f.readlines():
            if len(l.split('\n')[0].split(' ')) == 5:
                # x y z r
                lines.append([float(r) for r in l.split('\n')[0].split(' ')[:-1]])
            

    spheres = np.array(lines)*scale *scale_i

    return spheres


# 3. PARTICLE POSITION IN SPHERES
def getPartRelPos(obs_all, traj_part_):


    rel_pos        = np.full((traj_part_.shape[0],), False)   
    
    nb_parts_ini    = 10000
    ini_s           = 0 
    ini_e           = 0 
    while (ini_e < traj_part_.shape[0]):     
        ini_s = ini_e 
        ini_e = min(ini_s+nb_parts_ini, traj_part_.shape[0])
        curr_dist_part_sph   = ((np.expand_dims(traj_part_[ini_s:ini_e, :3], 1)- np.expand_dims( obs_all[:, :3], 0))**2).sum(axis=-1)
        rel_pos[ini_s:ini_e]= np.any((curr_dist_part_sph < ((obs_all[:, -1])**2)), axis=1)


    return rel_pos


# 4. Propagator

def update_W_elem(w2_elem, w4_elem, dd):

    for curr_w2_key in W2_COMB_IDX.keys():
        idx_ = W2_COMB_IDX[curr_w2_key]
        w2_elem[W2_IND_ELE[curr_w2_key]] += (dd[idx_[0]]*dd[idx_[1]]).sum(axis=0)

    for curr_comb_key in K_COMB_IDX:
        idx_ = K_COMB_IDX[curr_comb_key]
        w4_elem[K_IND_ELE[curr_comb_key]] += (dd[idx_[0]]*dd[idx_[1]]*dd[idx_[2]]*dd[idx_[3]]).sum(axis=0)

    return w2_elem, w4_elem 

def DT_from_elem(w2_elem, dt_vec): 
    DT = np.array([[w2_elem[0], w2_elem[1], w2_elem[2]], [w2_elem[1], w2_elem[3], w2_elem[4]], [w2_elem[2], w2_elem[4], w2_elem[5]]]) / 2 / dt_vec
    return DT

def KT_from_elem(w2_elem, w4_elem):

    K_elem  = np.zeros((15, w4_elem.shape[-1]))
    r       = 9 / (w2_elem[0]+w2_elem[3]+w2_elem[5])**2

    for curr_comb_key in K_COMB_IDX.keys(): 
        curr_w4         = w4_elem[K_IND_ELE[curr_comb_key]]
        curr_w4_idx     = K_COMB_IDX[curr_comb_key]
        curr_w2_key     = [int((curr_w4_idx[0]+1)*(curr_w4_idx[1]+1)), int((curr_w4_idx[2]+1)*(curr_w4_idx[3]+1)), int((curr_w4_idx[0]+1)*(curr_w4_idx[2]+1)), int((curr_w4_idx[1]+1)*(curr_w4_idx[3]+1)), int((curr_w4_idx[0]+1)*(curr_w4_idx[3]+1)), int((curr_w4_idx[1]+1)*(curr_w4_idx[2]+1))] 
        curr_w2_idx     = [W2_IND_ELE[curr_w2] for curr_w2 in curr_w2_key]

        Wijkl                                           = curr_w4 - w2_elem[curr_w2_idx[0]] * w2_elem[curr_w2_idx[1]] - w2_elem[curr_w2_idx[2]]*w2_elem[curr_w2_idx[3]] - w2_elem[curr_w2_idx[4]]*w2_elem[curr_w2_idx[5]]
        K_elem[K_IND_ELE[curr_comb_key]]      =  Wijkl * r 

    return K_elem

def ADC_from_DT(DT): 
    return np.trace(DT)/DT.shape[0]

def MK_from_KT(DT, K_elem):

    b_vals_dummy    = (np.random.rand(10))*1e9
    b_vecs_dummy    = np.random.rand(10, 3)
    b_vecs_dummy    /=  np.sqrt((b_vecs_dummy**2).sum(axis=1)).reshape(-1, 1)
    gtab_dummy      = gradient_table(b_vals_dummy, b_vecs_dummy, b0_threshold=b_vals_dummy.min()+1)
    dkimodel        = dki.DiffusionKurtosisModel(gtab_dummy)

    MK = np.zeros(DT.shape[-1])

    for i in range(DT.shape[-1]):
   
        eigval_DT, eigvec_DT    = np.linalg.eig(DT[:, :, i])

        idx                     = eigval_DT.argsort()[::-1]   
        eigval_DT               = eigval_DT[idx]
        eigvec_DT               = eigvec_DT[:,idx]

        dkiparams               = np.zeros(27)
        dkiparams[:3]           = eigval_DT
        dkiparams[3:12]         = eigvec_DT.flatten()
        
        for kkey in K_IND_ELE.keys():
            dkiparams[K_IND_ELE[kkey]+12] = K_elem[K_IND_ELE[kkey], i]

        dkifit_home = dki.DiffusionKurtosisFit(dkimodel, dkiparams)
        MK[i]       = dkifit_home.mk()

    return MK

# 5. SIGNAL
def calculate_DWI_signal_single_in_ex(G, idx_):

    ph = give_ph(G, idx_)
    DWI_real_in  = give_dwi(ph[init_pos])
    DWI_real_ex  = give_dwi(ph[np.invert(init_pos)])
    DWI_real_all = give_dwi(ph)
    return DWI_real_in, DWI_real_ex, DWI_real_all

def calculate_DWI_signal_mp_in_ex():

    DWI_real_in__   = np.zeros(nb_acq)
    DWI_real_ex__   = np.zeros_like(DWI_real_in__)
    DWI_real_all__  =  np.zeros_like(DWI_real_in__)

    with mp.Pool(nb_proc_max) as pool:
        mp_simu = pool.starmap(calculate_DWI_signal_single_in_ex, input_st )
        for ii, res_ in enumerate(mp_simu):
            DWI_real_in__[ii]   += res_[0]
            DWI_real_ex__[ii]   += res_[1]
            DWI_real_all__[ii]  += res_[2]           
            
    return DWI_real_in__, DWI_real_ex__, DWI_real_all__

def calculate_DWI_signal_single(G, idx_):
    ph = give_ph(G, idx_)
    DWI_real_all = give_dwi(ph)
    return DWI_real_all

def calculate_DWI_signal_mp():
    DWI_real_all__  = np.zeros(nb_acq)
    with mp.Pool(nb_proc_max) as pool:
        mp_simu = pool.starmap(calculate_DWI_signal_single, input_st )
        for ii, res_ in enumerate(mp_simu):
            DWI_real_all__[ii]  += res_
    return DWI_real_all__

def give_dwi(ph_sum):
    DWI_real  = np.sum(np.cos(ph_sum))     
    return DWI_real

def give_ph(G, idx_):

    ph = np.zeros(displ_x.shape[0])

    for part_ in range(displ_x.shape[0]):     
        
        x_p = displ_x[part_]
        y_p = displ_y[part_]
        z_p = displ_z[part_]

        ph[part_] =  np.sum((x_p[idx_]*G[0] +  y_p[idx_]*G[1] + z_p[idx_]*G[2]))%PI_2
            
    return ph

def get_displacement(traj):
    dx   = (traj[:, nb_inf::nb_inf] - traj[:, :1])
    dy   = (traj[:, nb_inf+1::nb_inf] - traj[:, 1:2])
    dz   = (traj[:, nb_inf+2::nb_inf] - traj[:, 2:3])
    return dx, dy, dz

def get_grad_vec(proto, duration, nb_step):
    if proto.type=="PGSE":
        Gvec, idx_vec   = get_grad_vec_PGSE(proto, duration, nb_step)
    elif proto.type=="STEAM":
        Gvec, idx_vec   = get_grad_vec_STEAM(proto, duration, nb_step)
    return Gvec, idx_vec

def get_grad_vec_PGSE(proto, duration, nb_step):

    t_vec   = np.linspace(0, duration, nb_step+1)
    dt      = t_vec[1] - t_vec[0]
    

    first_block_start_a   = ((proto.TE - proto.delta - proto.tdiff)/2).flatten()*1e-3
    first_block_end_a     = first_block_start_a + proto.delta.flatten()*1e-3
    second_block_start_a  = first_block_start_a + proto.tdiff.flatten()*1e-3
    second_block_end_a    = second_block_start_a + proto.delta.flatten()*1e-3

    Gdt_all = []
    idx_all = []

    nb_acq  = proto.G.size

    
    for i in range(nb_acq):
        first_block_start   = first_block_start_a[i]
        first_block_end     = first_block_end_a[i]
        second_block_start  = second_block_start_a[i]
        second_block_end    = second_block_end_a[i]

        # 1. Step in-between start and block 1
        idx_s_bl1   = np.where((t_vec[:-1] < first_block_start) & (t_vec[1:] > first_block_start))[0]
        
        # 2. Step in block 1
        idx_bl1     = np.where((t_vec[:-1] > first_block_start) & (t_vec[1:] < first_block_end))[0]
        
        # 3. Step in between block 1 and diffusion time
        idx_bl1_e   = np.where((t_vec[:-1] < first_block_end) & (t_vec[1:] > first_block_end))[0]
        
        #4. Step in between diffustion time and block 2
        idx_s_bl2   = np.where((t_vec[:-1] < second_block_start) & (t_vec[1:] > second_block_start))[0]
        
        #5. Step in block 2
        idx_bl2     = np.where((t_vec[:-1] > second_block_start) & (t_vec[1:] < second_block_end))[0]
        
        #6. Step in between block 2 and e)nd
        idx_bl2_e   = np.where((t_vec[:-1] < second_block_end) & (t_vec[1:] > second_block_end))[0]
        
        # Final
        dt_vec  = np.concatenate((t_vec[idx_s_bl1+1] - first_block_start ,np.zeros(idx_bl1.size) + dt, first_block_end - t_vec[idx_bl1_e], -(t_vec[idx_s_bl2+1] - second_block_start), np.zeros(idx_bl2.size) - dt, -(second_block_end - t_vec[idx_bl2_e])))
        idx     = np.concatenate((idx_s_bl1, idx_bl1, idx_bl1_e, idx_s_bl2, idx_bl2, idx_bl2_e))

        Gx      = dt_vec * proto.G[0,i]*proto.gdir[i][0]*GYRO*1e6
        Gy      = dt_vec * proto.G[0,i]*proto.gdir[i][1]*GYRO*1e6
        Gz      = dt_vec * proto.G[0,i]*proto.gdir[i][2]*GYRO*1e6
        
        Gdt_all.append((Gx, Gy, Gz))
        idx_all.append(idx)

    return Gdt_all, idx_all

def get_grad_vec_STEAM(proto, duration, nb_step):

    t_vec   = np.linspace(0, duration, nb_step+1)
    dt      = t_vec[1] - t_vec[0]

    nb_acq  = proto.G.size
    Gdt_all = []
    idx_all = []

    # Gradients 
    Gdiff   = (proto.gdir*proto.G.reshape(-1, 1))*GYRO*1e6
    Gc      = proto.cG*GYRO*1e6
    Gr      = proto.rG*GYRO*1e6

    # Diffusion gradient
    ts_td   = (proto.TE/2 - proto.delta_r - proto.delta_c - proto.tau1 - proto.delta).flatten()*1e-3
    te_td   = ts_td + proto.delta.flatten()*1e-3

    # Crusher gradient
    ts_tc   = te_td + proto.tau1.flatten()*1e-3
    te_tc   = ts_tc + proto.delta_c.flatten()*1e-3

    #Slice gradient 
    ts_tr   = te_tc 
    te_tr   = ts_tr + proto.delta_r.flatten()*1e-3

    # Slice gradient 2 
    ts_tr2  = te_tr + proto.tau_m.flatten() *1e-3
    te_tr2  = ts_tr2 + proto.delta_r.flatten()*1e-3

    # Crusher gradient 2
    ts_tc2  = te_tr2
    te_tc2  = ts_tc2 + proto.delta_c.flatten()*1e-3

    # Diffusion gradient 2
    ts_td2  = te_tc2 + proto.tau2.flatten() *1e-3
    te_td2  = ts_td2 + proto.delta.flatten() *1e-3

    for i in range(nb_acq):

        # Diffusion gradient
        idx_d1    = np.where((t_vec > ts_td[i]) & (t_vec < te_td[i]))[0]

        # Crusher gradient
        idx_c1    = np.where((t_vec > ts_tc[i]) & (t_vec < te_tc[i]))[0]

        #Slice gradient 
        idx_r1    = np.where((t_vec > ts_tr[i]) & (t_vec < te_tr[i]))[0]

        # Slice gradient 2 
        idx_r2    = np.where((t_vec > ts_tr2[i]) & (t_vec < te_tr2[i]))[0]

        # Crusher gradient 2[0]
        idx_c2    = np.where((t_vec > ts_tc2[i]) & (t_vec < te_tc2[i]))[0]

        # Diffusion gradient 2
        idx_d2    = np.where((t_vec > ts_td2[i]) & (t_vec < te_td2[i]))[0]

        idx = np.concatenate((idx_d1, idx_c1, idx_r1, idx_r2, idx_c2, idx_d2))

        Gx  = np.concatenate((np.full(idx_d1.size, dt*Gdiff[i, 0]), np.full(idx_c1.size, dt*Gc[i, 0]), np.full(idx_r1.size, dt*Gr[i, 0]), -np.full(idx_r2.size, dt*Gr[i, 0]), -np.full(idx_c2.size, dt*Gc[i, 0]), -np.full(idx_d2.size, dt*Gdiff[i, 0])))
        Gy  = np.concatenate((np.full(idx_d1.size, dt*Gdiff[i, 1]), np.full(idx_c1.size, dt*Gc[i, 1]), np.full(idx_r1.size, dt*Gr[i, 1]), -np.full(idx_r2.size, dt*Gr[i, 1]), -np.full(idx_c2.size, dt*Gc[i, 1]), -np.full(idx_d2.size, dt*Gdiff[i, 1])))
        Gz  = np.concatenate((np.full(idx_d1.size, dt*Gdiff[i, 2]), np.full(idx_c1.size, dt*Gc[i, 2]), np.full(idx_r1.size, dt*Gr[i, 2]), -np.full(idx_r2.size, dt*Gr[i, 2]), -np.full(idx_c2.size, dt*Gc[i, 2]), -np.full(idx_d2.size, dt*Gdiff[i, 2])))

        Gdt_all.append((Gx, Gy, Gz))
        idx_all.append(idx)

    return Gdt_all, idx_all



# 6. MAIN
if __name__=="__main__":

    """ 
    Example: 

    python PATH_TO_FOLDER_WITH_TRAJ PREFIX_OF_TRAJ_FILE --obstacle_file PATH_TO_SPHERE_LIST --scheme NAME_OF_SCHEME_FILE --propagator [0 or 1] 
    python dwi_signal_intra.py /home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/funnel/overlap_4/n1/ "N_50000_T_5000_" --scheme /home/localadmin/Documents/MCDC_perm_jas/Permeable_MCDS/results/funnel/overlap_4/n1/PGSE_21_dir_12_b.scheme --propagator 0
    
    
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('dir_folder', help = 'Folder of experiment. Folder to load trajectories and save results')
    parser.add_argument('exp_number', help = 'Experiment name. Prefix of trajectories file')
    parser.add_argument('--obstacle_file', help = 'Path to obstacle file', type=str, default="")
    parser.add_argument('--scheme', help = 'Scheme file MRI sequence', type=str, default="")
    parser.add_argument('--propagator', help = 'Flag for computing propagator',type=int, default=0)
    
    args = parser.parse_args()

    # 1. Parameters
    # A. General
    dir_folder      = args.dir_folder
    exp_number      = args.exp_number
    obstacle_file   = args.obstacle_file
    scheme_file     = args.scheme
    flag_in_ex      = False if obstacle_file=="" else True
    flag_signal     = False if scheme_file=="" else True
    flag_propagator = args.propagator

    print("Flags: In-ex = {} - Signal = {} -  Propagator = {} ".format(flag_in_ex, flag_signal, flag_propagator))

    dir_path    = dir_folder
    exp_path    = os.path.join(dir_path, exp_number)
    info        = read_info(os.path.join(dir_path, exp_number+'simulation_info.txt'))

    nb_proc_max     = 10
    
    nb_step         = info['T']
    nb_inf          = 3
    T_exp           = info["duration"]
    shape_          = (-1, (nb_step+1)*nb_inf)
    list_traj_files = [os.path.join(dir_path, a) for a in os.listdir(dir_path) if (('traj' in a ) and (exp_number in a))]
    nb_traj_file    = len(list_traj_files)

    nb_parts_all = 0 

    # B. Signal
    if flag_signal:

        path_scheme     = os.path.join(dir_path, scheme_file)
        proto           = Protocol(path_scheme, orig="text")
        nb_acq          = proto.tdiff.size
    
        save_path_r_all = os.path.join(exp_path+'_DWI.txt')
        DWI_real_all    = np.zeros(nb_acq)

        Gvec, idx_vec   = get_grad_vec(proto, T_exp, nb_step)
    
        input_st        = [(Gvec[j], idx_vec[j]) for j in range(nb_acq)]  

        if flag_in_ex:

            save_path_r_in  = os.path.join(exp_path+scheme_file.split('.')[0]+'_in_DWI.bfloat')
            save_path_r_ex  = os.path.join(exp_path+scheme_file.split('.')[0]+'_ex_DWI.bfloat')

            DWI_real_in     = np.zeros(nb_acq)
            DWI_real_ex     = np.zeros(nb_acq)


    # C. Propagator 
    if flag_propagator:
        idx_s       = 10
        N_prop      = 1000
        N_step_prop = 2000 

        time_step_simu      = np.arange(1, nb_step+1, dtype=int)* T_exp/nb_step
        time_step_prop_idx  = np.linspace(idx_s, nb_step-1, N_step_prop, dtype=int)

        nb_parts_per_file   = info["N"] / info["num_process"]
        nb_parts_up = 10000


        dist_mean           = np.sqrt(6*info["diffusivity"]*T_exp)*1e3 # mm
        dist_max            = (dist_mean*3)
        
        bins                   = np.linspace(0, dist_max, N_prop+1)

        propagator_all         = np.zeros((N_prop, N_step_prop))
            
        w2_elem_all            = np.zeros((6, nb_step))
        w4_elem_all            = np.zeros((15, nb_step))

        if flag_in_ex:

            propagator_in          = np.zeros((N_prop, N_step_prop))
            propagator_ex          = np.zeros((N_prop, N_step_prop))
            
            w2_elem_in             = np.zeros((6, nb_step))
            w4_elem_in             = np.zeros((15, nb_step))

            w2_elem_ex             = np.zeros((6, nb_step))
            w4_elem_ex             = np.zeros((15, nb_step))


    if flag_in_ex:
        nb_parts_in = 0 
        nb_parts_ex = 0 

        obs_list    = read_spheres(obstacle_file) *1e3

    # 2. Signal generation and propagator statistics
    for (curr_idx, curr_traj) in tqdm.tqdm(enumerate(list_traj_files), total=nb_traj_file):

        # traj_part_tmp = array.array('f')
        # with open(curr_traj, 'r') as f:
        #     lines = f.readlines()
        #     lines = [i for i in lines if i != '\n']
        #     traj_part_tmp = np.array(lines).astype(float)
        #     # traj_part_tmp.fromfile(f, os.path.getsize(curr_traj) // traj_part_tmp.itemsize)
        # traj_part   = np.reshape(np.array(traj_part_tmp), shape_)
        # print(traj_part.shape)

        traj_part_tmp = array.array('f')
        with open(curr_traj, 'rb') as f:
            traj_part_tmp.fromfile(f, os.path.getsize(curr_traj) // traj_part_tmp.itemsize)
        traj_part   = np.reshape(np.array(traj_part_tmp), shape_)
        
        del traj_part_tmp

        nb_parts_all += traj_part.shape[0]

        # Inital position
        if flag_in_ex:

            init_pos = getPartRelPos(obs_list, traj_part)
            
            nb_parts_in += init_pos.sum()
            nb_parts_ex += np.invert(init_pos).sum() 

        # Displacements
        displ_x, displ_y, displ_z  = get_displacement(traj_part)

        del traj_part

        # A. Signal
        if flag_signal:
            if flag_in_ex:

                DWI_real_in_, DWI_real_ex_, DWI_real_all_ = calculate_DWI_signal_mp_in_ex()
                DWI_real_in     += DWI_real_in_
                DWI_real_ex     += DWI_real_ex_
                DWI_real_all    += DWI_real_all_

                save_path_uni_all   = os.path.join(exp_path+scheme_file.split('.')[0]+'_'+str(curr_idx)+'_DWI.bfloat')
                save_path_uni_ex    = os.path.join(exp_path+scheme_file.split('.')[0]+'_'+str(curr_idx)+'_ex_DWI.bfloat')
                save_path_uni_in    = os.path.join(exp_path+scheme_file.split('.')[0]+'_'+str(curr_idx)+'_in_DWI.bfloat')

                # darray = array.array('f', DWI_real_in_)
                # with open(save_path_uni_in, 'wb') as f:
                #     darray.tofile(f)    
                    
                # darray = array.array('f', DWI_real_ex_)
                # with open(save_path_uni_ex, 'wb') as f:
                #     darray.tofile(f)    

                # darray = array.array('f', DWI_real_all_)
                # with open(save_path_uni_all, 'wb') as f:
                #     darray.tofile(f)  
   
            else:
                DWI_real_all_ = calculate_DWI_signal_mp()
                DWI_real_all    += DWI_real_all_

                # save_path_uni_all   = os.path.join(exp_path+'_'+str(curr_idx)+'_DWI.bfloat')
                # darray = array.array('f', DWI_real_all)
                # with open(save_path_uni_all, 'wb') as f:
                #     darray.tofile(f)                  
            
        # B. Propagator
        if flag_propagator:

            idx_s_up = 0
            idx_e_up = 0

            while (idx_e_up != (displ_x.shape[0] )):
            
                idx_s_up    = idx_e_up 
                idx_e_up    = np.min( (idx_s_up+nb_parts_up, displ_x.shape[0]))

                w2_elem_all, w4_elem_all  = update_W_elem(w2_elem_all, w4_elem_all, [displ_x[idx_s_up:idx_e_up], displ_y[idx_s_up:idx_e_up], displ_z[idx_s_up:idx_e_up]])

                if flag_in_ex:
                    idx_prop    = init_pos.copy()
                    idx_prop[:idx_s_up] = False
                    idx_prop[idx_e_up:] = False
                    w2_elem_in, w4_elem_in    = update_W_elem(w2_elem_in, w4_elem_in, [displ_x[idx_prop], displ_y[idx_prop], displ_z[idx_prop]])
                    
                    idx_prop            = np.invert(idx_prop)
                    idx_prop[:idx_s_up] = False
                    idx_prop[idx_e_up:] = False
                    w2_elem_ex, w4_elem_ex    = update_W_elem(w2_elem_ex, w4_elem_ex, [displ_x[idx_prop], displ_y[idx_prop], displ_z[idx_prop]])


            displ_ds    = np.sqrt(displ_x[:, time_step_prop_idx]**2 + displ_y[:, time_step_prop_idx]**2 + displ_z[:, time_step_prop_idx]**2)      

            for curr_step in range(displ_ds.shape[-1]):

                if flag_in_ex:
                    hist_in, _ = np.histogram(displ_ds[init_pos, curr_step], bins=bins)
                    propagator_in[:, curr_step] += hist_in

                    hist_ex, _ = np.histogram(displ_ds[np.invert(init_pos), curr_step], bins=bins)
                    propagator_ex[:, curr_step] += hist_ex

                    propagator_all[:, curr_step] += (hist_ex+hist_in)


                else:
                    hist_, _ = np.histogram(displ_ds[:, curr_step], bins=bins)
                    propagator_all[:, curr_step] += hist_

            del displ_ds

        del displ_x, displ_y, displ_z

    # Full signal
    if flag_signal:

        darray = array.array('f', DWI_real_all)
        with open(save_path_r_all, 'w') as f:
            #darray.tofile(f)  
            for DWI in DWI_real_all:
                f.write(f"{DWI}\n")
        if flag_in_ex:
            darray = array.array('f', DWI_real_in)
            with open(save_path_r_in, 'wb') as f:
                darray.tofile(f)    
                
            darray = array.array('f', DWI_real_ex)
            with open(save_path_r_ex, 'wb') as f:
                darray.tofile(f)    
        

    # Full prop
    if flag_propagator:
        
        if flag_in_ex:

            #intra
            w2_elem_in = w2_elem_in/ nb_parts_in
            w4_elem_in = w4_elem_in / nb_parts_in

            DT_prop_in = DT_from_elem(w2_elem_in, time_step_simu)
            KT_prop_in = KT_from_elem(w2_elem_in, w4_elem_in)
            ADC_prop_in= ADC_from_DT(DT_prop_in)
            MK_prop_in = MK_from_KT(DT_prop_in, KT_prop_in)
            
            np.save(os.path.join(exp_path+"propagator_DT_in.npy"), DT_prop_in)
            np.save(os.path.join(exp_path+"propagator_KT_in.npy"), KT_prop_in)
            np.save(os.path.join(exp_path+"propagator_ADC_in.npy"), ADC_prop_in)
            np.save(os.path.join(exp_path+"propagator_MK_in.npy"), MK_prop_in) 
            np.save(os.path.join(exp_path+"propagator_in.npy"), propagator_in)
            np.save(os.path.join(exp_path+"propagator_vec.npy"), bins)

            # Extra
            w2_elem_ex = w2_elem_ex / nb_parts_ex
            w4_elem_ex = w4_elem_ex / nb_parts_ex

            DT_prop_ex = DT_from_elem(w2_elem_ex, time_step_simu)
            KT_prop_ex = KT_from_elem(w2_elem_ex, w4_elem_ex)
            ADC_prop_ex= ADC_from_DT(DT_prop_ex)
            MK_prop_ex = MK_from_KT(DT_prop_ex, KT_prop_ex)
            
            np.save(os.path.join(exp_path+"propagator_DT_ex.npy"), DT_prop_ex)
            np.save(os.path.join(exp_path+"propagator_KT_ex.npy"), KT_prop_ex)
            np.save(os.path.join(exp_path+"propagator_ADC_ex.npy"), ADC_prop_ex)
            np.save(os.path.join(exp_path+"propagator_MK_ex.npy"), MK_prop_ex)    
            np.save(os.path.join(exp_path+"propagator_ex.npy"), propagator_ex)

        # All
        w2_elem_all = w2_elem_all / nb_parts_all
        w4_elem_all = w4_elem_all / nb_parts_all

        DT_prop_all = DT_from_elem(w2_elem_all, time_step_simu)
        KT_prop_all = KT_from_elem(w2_elem_all, w4_elem_all)
        ADC_prop_all= ADC_from_DT(DT_prop_all)
        MK_prop_all = MK_from_KT(DT_prop_all, KT_prop_all)
        
        np.save(os.path.join(exp_path+"propagator_DT.npy"), DT_prop_all)
        np.save(os.path.join(exp_path+"propagator_KT.npy"), KT_prop_all)
        np.save(os.path.join(exp_path+"propagator_ADC.npy"), ADC_prop_all)
        np.save(os.path.join(exp_path+"propagator_MK.npy"), MK_prop_all)    
        np.save(os.path.join(exp_path+"propagator.npy"), propagator_all)

    print('Done: File saved under {}'.format(exp_path))
    