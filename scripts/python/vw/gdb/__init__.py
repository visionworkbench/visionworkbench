#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__


from __future__ import absolute_import

import gdb
import re
import sys

from vw.gdb import types

DEBUG = False
def dprint(str):
    if DEBUG:
        print str

try:
    if not hasattr(sys, 'argv'):
        sys.argv = ['gdb']
    from IPython.Shell import IPShellEmbed
    _ipshell = IPShellEmbed(argv=[''],banner="Hello!",exit_msg="Goodbye!")
except ImportError:
    def noop():
        pass
    _ipshell = noop


pretty_printers_dict = None

def register_handler(types):
    for t in types:
        pretty_printers_dict[re.compile(t.PATTERN)] = t
        dprint('Registered %s for %s' % (t.__name__, t.PATTERN))

def build_dictionary():
    global pretty_printers_dict
    pretty_printers_dict = {}
    register_handler([getattr(types, i) for i in dir(types) if i[0] != '_' and i[0] == i[0].upper()])

def lookup(val):
    "Look-up and return a pretty-printer that can print val."

    if pretty_printers_dict is None:
        build_dictionary()

    # Get the type.
    typ = val.type

    # If it points to a reference, get the reference.
    if typ.code == gdb.TYPE_CODE_REF:
        typ = typ.target()

    # Get the unqualified type, stripped of typedefs.
    typ = typ.unqualified().strip_typedefs()

    # Get the type name.
    typename = typ.tag
    if typename == None:
        return None

    # Iterate over local dictionary of types to determine
    # if a printer is registered for that type.  Return an
    # instantiation of the printer if found.
    for function in pretty_printers_dict:
        m = function.match(typename)
        if m:
            handler = pretty_printers_dict[function]
            dprint('Matched %s for %s' % (handler.__name__, typename))
            if hasattr(handler, 'can_use'):
                if not handler.can_use(val):
                    return None
            return handler(typename, val)

    # Cannot find a pretty printer.
    return None

def register(obj):
    "Register pretty-printers with objfile obj."
    if obj == None:
        obj = gdb
    obj.pretty_printers.append(lookup)

__all__ = (register)
