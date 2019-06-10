#ifndef Y3_CLUSTER_DEFAULT_MODELS_HH
#define Y3_CLUSTER_DEFAULT_MODELS_HH


#include "models/a_cen_t.hh"
#include "models/a_mis_t.hh"
#include "models/del_sig_tom.hh"
#include "models/dv_do_dz_t.hh"
#include "models/hmb_t.hh"
#include "models/hmf_t.hh"
#include "models/lc_lt_t.hh"
#include "models/lo_lc_t.hh"
#include "models/mor_t2.hh"
#include "models/roffset_t.hh"
#include "models/sample_variance.hh"
#include "models/t_cen_t.hh"
#include "models/t_mis_t.hh"
#include "models/omega_z_sdss.hh"
#include "models/int_zo_zt_t.hh"
#include "models/pzsource_t.hh"
#include "models/pzsource_gaussian_t.hh"
#include "models/models.hh"

namespace y3_cluster {
    using DefaultModels = Models<MOR_t2,
                                 LO_LC_t,
                                 LC_LT_t,
                                 INT_ZO_ZT_t,
                                 ROFFSET_t,
                                 T_CEN_t,
                                 T_MIS_t,
                                 A_CEN_t,
                                 A_MIS_t,
                                 HMB_t,
                                 HMF_t,
                                 DEL_SIG_TOM,
                                 DV_DO_DZ_t,
                                 OMEGA_Z_SDSS>;
}

#endif
