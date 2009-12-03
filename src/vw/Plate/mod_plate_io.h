#ifndef __VW_PLATE_MOD_PLATE_IO_H__
#define __VW_PLATE_MOD_PLATE_IO_H__

#include <httpd.h>

#ifdef __cplusplus
extern "C" {
#endif

int  mod_plate_handler(request_rec *r);
int  mod_plate_status(request_rec *r, int flags);
void mod_plate_child_init(apr_pool_t *pchild, server_rec *s);

#ifdef __cplusplus
}
#endif

#endif
