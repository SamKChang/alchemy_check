
# coding: utf-8

# In[ ]:

import qctoolkit as qtk
import numpy as np
import os
import glob


# In[ ]:

# construct dummy base molecule template
base = qtk.Molecule()
base.celldm = [20, 10, 10, 0, 0, 0]
base.isolated = True
hf = base.copy()
hf.build([[1, 5., 5., 5.], [9, 6., 5., 5.]])
hf.name = 'hf_ref'
hcl = base.copy()
hcl.build([[1, 5., 5., 5.], [17, 6., 5., 5.]])
hcl.name = 'hcl_tar'

# basic setup
prefix = 'production_shifted_'
ref_root = os.path.abspath('%srefs' % prefix)

qmsettings = {
    'program': 'cpmd',
    'omp': 2,
    'cutoff': 200,
    'wf_convergence': 1E-7,
    #'threads': 4
}


# In[ ]:

# construct QM jobs
inp_refs = []
inp_prds = []
for d in np.arange(0.5, 2.4, 0.1):
    # construct HF reference run with RESTART file
    s = d - 1.
    d_str = '%02d' % np.round(d * 10)
    tmp = hf.copy()
    tmp.stretch(1, [0, 1], s)
    tmp.name = 'hf_ref%s' % d_str
    inp = qtk.QMInp(tmp, save_restart=True, **qmsettings)
    inp_refs.append(inp)

    # construct HCl target run
    tmp2 = hcl.copy()
    tmp2.stretch(1, [0, 1], s)
    tmp2.name = 'hcl_ref%s' % d_str
    inp = qtk.QMInp(tmp2, **qmsettings)
    inp_refs.append(inp)
    
    # construct HCl prediction run with alchemical prediction setting
    prd =  hcl.copy()
    prd.stretch(1, [0, 1], s)
    prd.name = 'hcl_prd%s' % d_str
    rst_file = '%s/hf_ref%s/RESTART.1' % (ref_root, d_str)
    inp = qtk.QMInp(
        prd, 
        scf_step=1,
        save_restart=True,
        restart=True,
        restart_wavefunction_file=rst_file,
        **qmsettings)
    inp_prds.append(inp)

# run all reference/target
qtk.qmRunAll(inp_refs, ref_root)
# run all prediction
qtk.qmRunAll(inp_prds, '%sprds' % prefix)


# In[ ]:

# load computed results

refs = []
for f in sorted(glob.glob('%srefs/hf*/*.out' % prefix)):
    refs.append(qtk.QMOut(f, program='cpmd'))
tars = []
for f in sorted(glob.glob('%srefs/hcl*/*.out' % prefix)):
    tars.append(qtk.QMOut(f, program='cpmd'))
prds = []
for f in sorted(glob.glob('%sprds/hcl*/*.out' % prefix)):
    prds.append(qtk.QMOut(f, program='cpmd'))

# extract bond length of each system
R = np.array([o.molecule.R[1,0] for o in refs])
    
# extract total energy of each system
re = np.array([o.Et for o in refs])
te = np.array([o.Et for o in tars])
pe = np.array([o.Et for o in prds])

# shift minimum to zero
res = re - re.min()
tes = te - te.min()
pes = pe - pe.min()

qtk.save(
    [
        [refs, tars, prds],
        R,
        [re, te, pe],
        [res, tes, pes]
    ],
    '%sresults.pkl' % prefix
)

