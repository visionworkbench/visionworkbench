#ifndef __VW_PLATE_MOD_PLATE_IO_H__
#define __VW_PLATE_MOD_PLATE_IO_H__

#include <httpd.h>

#ifdef __cplusplus
extern "C" {
#endif

int  mod_plate_handler(request_rec *r);
int  mod_plate_status(request_rec *r, int flags);
void mod_plate_child_init(apr_pool_t *pchild, server_rec *s);

typedef struct {
  const char *name;
  int level; // actually a vw::MessageLevel, but we can't include Log.h here.
} rule_entry;

typedef struct {
  const char *rabbit_ip;
  apr_array_header_t *rules; // This holds rule_entries
} plate_config;

const plate_config* get_plate_config(const server_rec* s);
plate_config* get_plate_config_mutable(server_rec* s);

#ifdef __cplusplus
}
#endif

#endif
