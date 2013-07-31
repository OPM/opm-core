/* this file is written by the build system; don't include it anywhere else */
#include <project-version.h>

/* declaration of the externally visible symbol */
#include <opm/core/version.h>

/* initialization */
const char* const opm_core_version = PROJECT_VERSION;
