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


// Types of template decl:
//
// [class/struct]
//     template <class PixelT> struct ValueEdgeExtension;
//   becomes
//     template struct ValueEdgeExtension<PixelT>;
//
// [free function]
//     template <class ImageT>
//     EdgeExtensionView<ImageT, ConstantEdgeExtension> edge_extend(const ImageViewBase<ImageT>&);
//   becomes
//     template
//     EdgeExtensionView<ImageT, ConstantEdgeExtension> edge_extend(const ImageViewBase<ImageT>&);

if (!String.prototype.pad)
{
    String.prototype.pad = function(l, s, t){
        return s || (s = " "), (l -= this.length) > 0 ? (s = new Array(Math.ceil(l / s.length)
            + 1).join(s)).substr(0, t = !t ? l : t == 1 ? 0 : Math.ceil(l / 2))
            + this + s.substr(0, l - t) : this;
    };
}

function removeSubstring(s, t) {
    i = s.indexOf(t);
    r = "";
    if (i == -1) return s;
    r += s.substring(0,i) + removeSubstring(s.substring(i + t.length), t);
    return r;
}

if (!Array.prototype.map)
{
  Array.prototype.map = function(fun /*, thisp*/)
  {
    var len = this.length;
    if (typeof fun != "function")
      throw new TypeError();

    var res = new Array(len);
    var thisp = arguments[1];
    for (var i = 0; i < len; i++)
    {
      if (i in this)
        res[i] = fun.call(thisp, this[i], i, this);
    }

    return res;
  };
}



function ModuleList()
{
    this.modules     = {};

    this.RECORD      = "Record";
    this.FUNC_FREE   = "Free";
    this.FUNC_MEMBER = "Member";

    this.module_name = function(d) {
        var file = d.loc.file;
        var idx  = d.loc.file.indexOf("src/vw/");
        if (idx == -1)
            throw Error("Could not identify module of " + file);
        idx += 7;

        return d.loc.file.substring(idx, file.indexOf("/", idx));
    };

    this.parseDecl = function(d) {
        if ((obj = this.Func(d)) == null && (obj = this.Record(d)) == null)
            return null;

        obj.loc = "// " + d.loc.file + ":" + d.loc.line;
        return obj;
    }

    this.Record = function(d) {
        if (!d.type || !d.type.name || !(d.type.kind == "struct" || d.type.kind == "class"))
            return null;

        if (d.type.template && d.type.template.arguments)
        {
            var args = d.type.template.arguments;
            for (var i in args)
            {
                if (args[i].isTypename) {
                    return {
                        name: d.type.name,
                        decl: "template " + d.type.kind.pad(7, " ", 1) + d.type.name + ";",
                        type: this.RECORD
                    }
                }
            }
        }
        return null;
    };

    this.Func = function(d) {
        if (!d.isFunction || !d.template || d.memberOf || !d.type.type.name)
            return null;

        var parts = d.name.split("(")
        parts[0] += " < " + d.template.map(function(v,i,a) {return v.name}).join() + " > ";
        var newname = removeSubstring(parts.join("("), "typename ");
        var decl = "template " + removeSubstring(d.type.type.name, "typename ") + " " + newname + ";";
        return {
            name: newname,
            decl: decl,
            type: this.FUNC_FREE
        };
    };

    return true;
}


ModuleList.prototype.add_decl = function(d) {
    var obj = this.parseDecl(d);

    if (!obj)
        return;

    var name = this.module_name(d);
    if (this.modules[name] == null)
        this.modules[name] = {};
    if (this.modules[name][obj.type] == null)
        this.modules[name][obj.type] = [];
    this.modules[name][obj.type].push(obj);
}

var list = new ModuleList();

function process_decl(d)
{
    if (!d.name || d.name.substring(0, 4) != "vw::")
        return;

    //if (d.name.indexOf("compound_select_channel") == -1)
    //    return;
    //print(d);
    //return;

    //if (d.name.indexOf("ArgValDifferenceFunctor::operator") == -1)
    //    return;

    //print(d);
    //if (d.name.indexOf("edge_extend") == -1)
    //    return;

    list.add_decl(d);
}

function sort_module_list(a,b) {
    return a.type < b.type ? -1 :
           a.type > b.type ?  1 :
           a.name < b.name ? -1 :
           a.name > b.name ?  1 :
           0;
}

function sort_name(a,b) {
    return a.name < b.name ? -1 :
           a.name == b.name ? 0 :
           1;
}

// try to back up a file before overwriting it. don't rely on this; it's fragile.
function write_file_bak(file, data)
{
    if (options['no-backup']) {
        write_file(file, data);
        return;
    }
    var old_data;
    try {old_data = read_file(file);}
    catch (e) {
        if (e.name != "Error" || e.message.substring(0,29) != "read_file: error opening file")
            throw e;
    }

    if (old_data)
    {
        if (old_data != data) {
          print("File " + file + " changed. saving backup.");
          write_file(file + "~", old_data);
          write_file(file, data);
        }
    }
    else
      write_file(file, data);
}

if (!options['output-path'])
    throw Error("please set OUTPUT_PATH");

function input_end()
{
    var fn;
    var excl;

    for (var module_name in list.modules)
    {
        var module = list.modules[module_name];

        excl = "";
        try {excl = read_file(options['output-path'] + "/" + module_name + "/tests/TestInstantiateExclusion.txt")}
        catch (e) {/* swallow. assume there aren't any exclusions. shame on me. */ }

        excl = removeSubstring(excl, " ");

        for (var type_name in module) {
            var type = module[type_name];
            type.sort(sort_name);
            var data = "";
            for (var d in type) {
                var d = type[d];
                //data += d.loc + "\n" + d.decl + "\n";
                if (excl.indexOf(removeSubstring(d.decl, " ")) != -1) // wow. linear search, huh? lame.
                    data += "//"
                data += d.decl + "\n";
            }
            if (data.length > 0)
              write_file_bak(options['output-path'] + "/" + module_name + "/tests/TestInstantiate" + type_name + "List.hh", data);
        }
    }
}
