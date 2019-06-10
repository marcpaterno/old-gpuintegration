//in normal circumstance, you don't need to modify this file
#ifndef Y3_CLUSTER_MODELS_HH
#define Y3_CLUSTER_MODELS_HH

namespace y3_cluster
{
    template<typename MOR_,
             typename LO_LC_,
             typename LC_LT_,
             typename ZO_ZT_,
             typename ROFFSET_,
             typename T_CEN_,
             typename T_MIS_,
             typename A_CEN_,
             typename A_MIS_,
             typename HMB_,
             typename HMF_,
             typename DEL_SIG_,
             typename DV_DO_DZ_,
             typename OMEGA_Z_>
    struct Models {
        using MOR = MOR_;
        using LO_LC = LO_LC_;
        using LC_LT = LC_LT_;
        using ZO_ZT = ZO_ZT_;
        using ROFFSET = ROFFSET_;
        using T_CEN = T_CEN_;
        using T_MIS = T_MIS_;
        using A_CEN = A_CEN_;
        using A_MIS = A_MIS_;
        // Halo Mass Bias and Halo Mass Function
        using HMB = HMB_;
        using HMF = HMF_;
        using DEL_SIG = DEL_SIG_;
        using DV_DO_DZ = DV_DO_DZ_;
        using OMEGA_Z = OMEGA_Z_;
    };

}

#endif
