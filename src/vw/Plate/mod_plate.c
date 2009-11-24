#include <httpd.h>
#include <ap_mpm.h>
#include <http_config.h>
#include <mod_status.h>

/* Forward declarations (You got c++ in my c!) */
#include "mod_plate_io.h"

static void plate_register_hooks(apr_pool_t *p)
{
  int threaded;

  if (ap_mpm_query(AP_MPMQ_IS_THREADED, &threaded) != APR_SUCCESS) {
    fprintf(stderr, "mod_plate could not query mpm: we have a problem. Refusing to register!");
    return;
  }

  if (threaded) {
    fprintf(stderr, "Refusing to start mod_plate inside a threaded MPM. We don't know how our threads will interact with apache's.");
    return;
  }

  mod_plate_init();
  ap_hook_handler(mod_plate_handler, NULL, NULL, APR_HOOK_MIDDLE);
  APR_OPTIONAL_HOOK(ap, status_hook, mod_plate_status, NULL, NULL, APR_HOOK_MIDDLE);
}

/* Dispatch list for API hooks */
module AP_MODULE_DECLARE_DATA plate_module = {
  STANDARD20_MODULE_STUFF, 
  NULL,                  /* create per-dir    config structures */
  NULL,                  /* merge  per-dir    config structures */
  NULL,                  /* create per-server config structures */
  NULL,                  /* merge  per-server config structures */
  NULL,                  /* table of config file commands       */
  plate_register_hooks  /* register hooks                      */
};

