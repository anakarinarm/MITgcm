C CPP options file for the mom_common package.

#ifndef MOM_COMMON_OPTIONS_H
#define MOM_COMMON_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_MOM_COMMON
#define COSINEMETH_III
#undef ISOTROPIC_COS_SCALING
#undef ALLOW_LEITH_QG
#undef ALLOW_SMAG_3D
#undef ALLOW_3D_VISCAH
#undef ALLOW_3D_VISCA4
#undef ALLOW_BOTTOMDRAG_ROUGHNESS
#undef MOM_USE_OLD_DEEP_VERT_ADV

C Separate explicit bottom-drag tendencies from aggregate dissipation.
C This is compatible with the configured selectImplicitDrag=0.
#define ALLOW_MOM_TEND_EXTRA_DIAGS
#endif

#endif
