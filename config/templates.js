// Types of template decl:
//
// [class/struct]
//    template <class PixelT> struct ValueEdgeExtension;
//    becomes
//    template struct ValueEdgeExtension<PixelT>;

function hasTmplArg(d)
{
    if (d.type.template && d.type.template.arguments)
    {
        var args = d.type.template.arguments;
        for (var i in args)
        {
            if (args[i].isTypename)
                return true;
        }
    }
    return false;
}

String.prototype.pad = function(l, s, t){
    return s || (s = " "), (l -= this.length) > 0 ? (s = new Array(Math.ceil(l / s.length)
        + 1).join(s)).substr(0, t = !t ? l : t == 1 ? 0 : Math.ceil(l / 2))
        + this + s.substr(0, l - t) : this;
};

var records = []
function process_decl(d)
{
    if (!d.name || d.name.substring(0, 4) != "vw::")
        return;

    if (!d.type || !d.type.name || !(d.type.kind == "struct" || d.type.kind == "class"))
        return;

    if (!hasTmplArg(d))
        return;

    records.push({
        loc:  "// " + d.loc.file + ":" + d.loc.line,
        module:  d.loc.file.substring(16, d.loc.file.indexOf("/", 16)),
        name: d.type.name,
        decl: "template " + d.type.kind.pad(7, " ", 1) + d.type.name + ";"
    });

}

function input_end()
{
    records.sort(function(a,b) {
            return a.module < b.module ? -1 :
                   a.module > b.module ?  1 :
                   a.name   < b.name   ? -1 :
                   a.name   > b.name})

    var last = "";
    var data = "";
    for (var d in records)
    {
        var d = records[d];
        if (d.module != last)
        {
            print("Switching to " + d.module + " from " + last);
            if (last)
                write_file("Instantiate" + last + ".h", data);
            last = d.module;
            data = "";
        }
        //data += d.loc + "\n" + d.decl + "\n";
        data += d.decl + "\n";
    }
    write_file("Instantiate" + last + ".h", data);
}
