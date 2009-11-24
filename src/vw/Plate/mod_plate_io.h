#ifndef __VW_PLATE_MOD_PLATE_IO_H__
#define __VW_PLATE_MOD_PLATE_IO_H__

#include <httpd.h>

#ifdef __cplusplus
extern "C" {
#endif

void mod_plate_init();
void mod_plate_destroy();
int  mod_plate_handler(request_rec *r);
int  mod_plate_status(request_rec *r, int flags);

#ifdef __cplusplus
}
#endif

#endif
