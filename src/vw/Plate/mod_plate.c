// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <httpd.h>
#include <http_config.h>
#include <mod_status.h>

#include <vw/Plate/common.h>

/* Forward declarations (You got c++ in my c!) */
#include "mod_plate.h"

module AP_MODULE_DECLARE_DATA plate_module;

static void plate_register_hooks(apr_pool_t *p)
{
  // warning go away
  (void)p;

  ap_hook_handler(mod_plate_handler, NULL, NULL, APR_HOOK_MIDDLE);
  ap_hook_child_init(mod_plate_child_init, NULL, NULL, APR_HOOK_MIDDLE);
  ap_hook_post_config(mod_plate_post_config, NULL, NULL, APR_HOOK_MIDDLE);
  APR_OPTIONAL_HOOK(ap, status_hook, mod_plate_status, NULL, NULL, APR_HOOK_MIDDLE);
}

static const char* handle_rule(cmd_parms* cmd, void* cfg, const char* raw_level, const char* name) {
  // if an exception is thrown from inside this function (anywhere) it will crash at throw-time.

  // warning go away
  (void)cfg;

  int level;

  if (strcmp(raw_level, "ErrorMessage") == 0)
    level = 0;
  else if (strcmp(raw_level, "WarningMessage") == 0)
    level = 10;
  else if (strcmp(raw_level, "InfoMessage") == 0)
    level = 20;
  else if (strcmp(raw_level, "DebugMessage") == 0)
    level = 30;
  else if (strcmp(raw_level, "VerboseDebugMessage") == 0)
    level = 40;
  else if (strcmp(raw_level, "EveryMessage") == 0)
    level = 100;
  else if (strcmp(raw_level, "*") == 0)
    level = 100;
  else
    return "Illegal log level";

  plate_config *conf = get_plate_config_mutable(cmd->server);
  rule_entry *rule   = (rule_entry*)apr_array_push(conf->rules);

  if (!name)
    name = "plate.apache";

  rule->name  = name;
  rule->level = level;

  return NULL;
}

static const char* handle_alias(cmd_parms* cmd, void* cfg, const char* name, const char* id) {
  plate_config *conf = get_plate_config_mutable(cmd->server);

  // warning go away
  (void)cfg;

  intptr_t id_i = atoi(id);
  if (id_i == 0)
    return "Illegal platefile id: 0";

  // XXX: intptr_t guarantees that a pointer can be cast to it and back without
  // changing it, but it doesn't guarantee that an intptr_t can be cast to a
  // pointer and back without changing it. So this code is technically
  // undefined. Sigh.
  apr_table_setn(conf->alias, name, (const char*)id_i);

  return NULL;
}

#define ADD_STRING_CONFIG(key, check) \
  static const char* handle_ ## key(cmd_parms* cmd, void* null, const char* arg) {\
    /* warning go away*/ \
    (void)null; \
    plate_config *cfg = get_plate_config_mutable(cmd->server);\
    cfg->key = arg;\
    const char *error = check(cfg->key);\
    if (error)\
      return error;\
    ap_set_module_config(cmd->server->module_config, &plate_module, (void*)cfg);\
    return NULL;\
  }

#define ADD_INT_CONFIG(key, check) \
  static const char* handle_ ## key(cmd_parms* cmd, void* null, const char* arg) {\
    /* warning go away*/ \
    (void)null; \
    plate_config *cfg = get_plate_config_mutable(cmd->server);\
    cfg->key = atoi(arg);\
    const char *error = check(cfg->key);\
    if (error)\
      return error;\
    ap_set_module_config(cmd->server->module_config, &plate_module, (void*)cfg);\
    return NULL;\
  }

#define ADD_FLAG_CONFIG(key)\
  static const char* handle_ ## key(cmd_parms* cmd, void* null, int on) {\
    /* warning go away*/ \
    (void)null; \
    plate_config *cfg = get_plate_config_mutable(cmd->server);\
    cfg->key = on ? 1 : 0;\
    ap_set_module_config(cmd->server->module_config, &plate_module, (void*)cfg);\
    return NULL;\
  }

const char* is_ip_address(const char* arg) {
    // really stupid heuristic
    if (!arg)
        return "Expected an IP address: cannot be null";
    if (arg[0] < '0' || arg[0] > '9')
        return "Expected an IP address: should have a leading number";
    return NULL;
}

const char* is_bare_exchange(const char* arg) {
    if (!arg)
        return "Expected a bare exchange name: cannot be null";
    if (strchr(arg, '.'))
        return "Expected a bare exchange name: periods delimit namespaces!";
    return NULL;
}

const char* is_platefile_id(const char* arg) {
  if (!arg)
    return "Expected a platefile id: cannot be null";

  char *end;
  if (strtol(arg, &end, 10) == 0)
    return "Expected a platefile id: got something else?";
  if (*end != '\0')
    return "Expected a platefile id: got a number and some junk";
  return NULL;
}

const char* is_servername(const char* arg) {
  if (!arg)
    return "Expected a servername: cannot be null";
  if (strncmp(arg, "http", 4) != 0)
    return "Expected a servername: it should probably begin with http";
  return NULL;
}

const char* is_gt_zero(int arg) {
  if (arg < 1)
    return "Expected a number greater than zero";
  return NULL;
}

ADD_STRING_CONFIG(rabbit_ip, is_ip_address);
ADD_STRING_CONFIG(index_exchange, is_bare_exchange);
ADD_STRING_CONFIG(dem_id, is_platefile_id);
ADD_STRING_CONFIG(servername, is_servername);
ADD_INT_CONFIG(index_timeout, is_gt_zero);
ADD_INT_CONFIG(index_tries, is_gt_zero);
ADD_FLAG_CONFIG(unknown_resync);
ADD_FLAG_CONFIG(use_blob_cache);

static const command_rec my_cmds[] = {
  AP_INIT_TAKE1("PlateRabbitMQIP",    handle_rabbit_ip,      NULL, RSRC_CONF, "The IP of the rabbitmq server"),
  AP_INIT_TAKE1("PlateIndexExchange", handle_index_exchange, NULL, RSRC_CONF, "The rabbitmq exchange on which to look for the index server"),
  AP_INIT_TAKE1("PlateDemID",         handle_dem_id,         NULL, RSRC_CONF, "The platefile ID of the DEM layer to use"),
  AP_INIT_TAKE1("PlateServerName",    handle_servername,     NULL, RSRC_CONF, "The servername to use for the server in the WTML"),
  AP_INIT_TAKE12("PlateLogRule",      handle_rule,           NULL, RSRC_CONF, "A log rule to add to the vw::RuleSet"),
  AP_INIT_TAKE1("PlateIndexTimeout",  handle_index_timeout,  NULL, RSRC_CONF, "How long to wait for the index_server to respond"),
  AP_INIT_TAKE1("PlateIndexTries",    handle_index_tries,    NULL, RSRC_CONF, "How many times to try talking to the index_server"),
  AP_INIT_TAKE2("PlateAlias",         handle_alias,          NULL, RSRC_CONF, "Name-to-platefile_id mappings"),
  AP_INIT_FLAG("PlateUnknownResync",  handle_unknown_resync, NULL, RSRC_CONF, "Should we resync the platefile list when someone asks for an unknown one?"),
  AP_INIT_FLAG("PlateBlobCache",      handle_use_blob_cache, NULL, RSRC_CONF, "Should the blob cache be used?"),
  { NULL }
};

static void* create_plate_config(apr_pool_t* p, server_rec* s) {
  // warning go away
  (void)s;

  plate_config* conf = (plate_config*)apr_pcalloc(p, sizeof(plate_config));
  conf->rabbit_ip      = "127.0.0.1";
  conf->index_exchange = DEV_INDEX_BARE;
  conf->dem_id         = NULL;
  conf->servername     = NULL;
  conf->index_timeout  = 3000;
  conf->index_tries    = 3;
  conf->rules = apr_array_make(p, 2, sizeof(rule_entry));
  conf->alias = apr_table_make(p, 4);
  conf->unknown_resync = 1;
  conf->use_blob_cache = 1;
  return conf;
  // This is the default config file
#if 0
  PlateRabbitMQIP 127.0.0.1
  PlateIndexExchange index
  PlateLogRule DebugMessage plate.apache
  PlateIndexTimeout 3000
  PlateIndexTries 3
  PlateUnknownResync on
  PlateBlobCache on
#endif

// these keys not set by default, but here are examples of possible valid ones
#if 0
  PlateDemID 123456789
  PlateServerName http://198.10.124.50
  PlateAlias hirise 123456789
#endif
}

plate_config* get_plate_config_mutable(server_rec* s) {
  return (plate_config*)ap_get_module_config(s->module_config, &plate_module);
}

const plate_config* get_plate_config(const server_rec* s) {
  return (const plate_config*)ap_get_module_config(s->module_config, &plate_module);
}

/* Dispatch list for API hooks */
module AP_MODULE_DECLARE_DATA plate_module = {
  STANDARD20_MODULE_STUFF,
  NULL,                  /* create per-dir    config structures */
  NULL,                  /* merge  per-dir    config structures */
  create_plate_config,   /* create per-server config structures */
  NULL,                  /* merge  per-server config structures */
  my_cmds,               /* table of config file commands       */
  plate_register_hooks  /* register hooks                      */
};
