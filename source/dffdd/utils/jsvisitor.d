module dffdd.utils.jsvisitor;

import std.json;
import std.path;
import std.file;
import std.exception;
import std.stdio;

struct VisitorData
{
    JSONValue* _root;
    JSONValue* _entry;
    bool isDir;
    string dir;
    string[] parentKeys;
    uint oversampling;
    ptrdiff_t offsetTX;
    real samplingFreq;
    JSONValue[string] fields;

    ref JSONValue root() @property { return *_root; }
    ref JSONValue entry() @property { return *_entry; }
    string currKey() @property { return parentKeys[$-1]; }
    string txPrefix() @property { return root["send_fn_prefix"].str; }
    string rxPrefix() @property { return root["recv_fn_prefix"].str; }
    string ext() @property { return root["data_fn_ext"].str; }
    string txFileName() @property { return buildPath(dir, txPrefix ~ currKey ~ ext); }
    string rxFileName() @property { return buildPath(dir, rxPrefix ~ currKey ~ ext); }
}


enum bool isVisitor(T) = is(typeof((T v){
    JSONValue e;
    VisitorData d;

    v = v.onDirectory(d, e);
    v.onFile(d, e);
}));


void visitJSON(V)(ref JSONValue root, V visitor)
{
    VisitorData d;
    d._root = &root;
    d._entry = &root;
    d.dir = d.root["prefix"].str;
    d.parentKeys ~= root["name"].str;

    visitJSONImpl(d, visitor);
}


void visitJSONImpl(V)(VisitorData data, V visitor)
{
    enforce(data.entry.type == JSON_TYPE.OBJECT);

    foreach(k, ref v; data.entry.object)
    {
        if(v.type != JSON_TYPE.OBJECT) continue;
        if("ignore" in v.object) continue;

        if("is_dir" in v.object){
            VisitorData d = data;
            d._entry = &v;
            d.isDir = true;
            d.dir = buildPath(d.dir, k);
            d.parentKeys ~= k;

            visitJSONImpl(d, visitor.onDirectory(d));
        }else{
            VisitorData d = data;
            d._entry = &v;
            d.isDir = false;
            d.parentKeys ~= k;
            d.oversampling = cast(uint)v["oversampling"].integer;
            d.samplingFreq = v["fs"].floating;
            d.offsetTX = v["offset"].integer;
            d.fields = v.object;

            visitor.onFile(d);
        }
    }
}
