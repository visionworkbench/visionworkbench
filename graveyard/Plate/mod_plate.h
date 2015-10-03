// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_PLATE_MOD_PLATE_H__
#define __VW_PLATE_MOD_PLATE_H__

// forward decls
struct server_rec;
struct request_rec;
struct apr_pool_t;
struct apr_array_header_t;
struct apr_table_t;

#ifdef __cplusplus
extern "C" {
#endif

int  mod_plate_handler(request_rec *r);
int  mod_plate_status(request_rec *r, int flags);
void mod_plate_child_init(apr_pool_t *pchild, server_rec *s);
int  mod_plate_post_config(apr_pool_t *pconf, apr_pool_t *plog, apr_pool_t *ptemp, server_rec *s);

typedef struct {
  const char *name;
  int level; // actually a vw::MessageLevel, but we can't include Log.h here.
} rule_entry;

typedef struct {
  const char *index_url;
  const char *dem_id;
  const char *servername;
  int index_timeout;
  int index_tries;
  int unknown_resync;
  int use_blob_cache;
  apr_array_header_t *rules; // This holds rule_entries
  apr_table_t *alias;        // key is name, value is id, an int stored as a const char*
} plate_config;

const plate_config* get_plate_config(const server_rec* s);
plate_config* get_plate_config_mutable(server_rec* s);

#ifdef __cplusplus
}
#endif

#endif
