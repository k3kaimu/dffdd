/**
任意の構造体やクラスをJSONに変換したり、その逆変換を行う
*/

module dffdd.utils.json;

import std.algorithm;
import std.array;
import std.complex;
import std.conv;
import std.exception;
import std.json;
import std.process;
import std.range;
import std.string;
import std.traits;
import std.typecons;
import std.typetuple;
import std.variant;
import std.functional;



private template _StaticIota(size_t M, size_t N)
if(M <= N)
{
    static if(M == N)
        alias _StaticIota = TypeTuple!();
    else
        alias _StaticIota = TypeTuple!(M, _StaticIota!(M+1, N));
}


class JSONException : Exception
{
    this(string msg, string file = null, size_t line = 0)
    {
        super(msg, file, line);
    }
}


/**

*/
template JSONEnv(alias overloads)
{
    /**

    */
    void fromJSONValue(T)(JSONValue json, ref T dst)
    {
        static if(is(typeof(overloads.fromJSONValueImpl(json, dst))))
            overloads.fromJSONValueImpl(json, dst);
        else static if(is(typeof(T.fromJSONValueImpl(json)) : T))
            dst = T.fromJSONValueImpl(json);
        else
            fromJSONValueImpl(json, dst);
    }


    ///
    JSONValue toJSONValue(T)(auto ref T t)
    {
        static if(is(typeof(overloads.toJSONValueImpl(t)) == JSONValue))
            return overloads.toJSONValueImpl(forward!t);
        else static if(is(typeof(t.toJSONValueImpl()) == JSONValue))
            return t.toJSONValueImpl();
        else
            return toJSONValueImpl(forward!t);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(is(T == JSONValue))
    {
        return value;
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(is(T == typeof(null)))
    out(result){
        assert(result.type == JSONType.null_);
    }
    do{
        return JSONValue(null);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(isSomeString!T)
    out(result){
        assert(result.type == JSONType.string);
    }
    do{
        return JSONValue(value.to!string);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(isUnsigned!T && isIntegral!T)
    out(result){
        assert(result.type == JSONType.uinteger);
    }
    do{
        return JSONValue(value);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(isSigned!T && isIntegral!T)
    out(result){
        assert(result.type == JSONType.integer);
    }
    do{
        return JSONValue(value);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(is(T == bool))
    out(result){
        assert(result.type == JSONType.true_ || result.type == JSONType.false_);
    }
    do{
        return JSONValue(value);
    }


    ///
    JSONValue toJSONValueImpl(T)(T value)
    if(isFloatingPoint!T)
    out(result){
        assert(result.type == JSONType.float_);
    }
    do{
        return JSONValue(value);
    }


    ///
    JSONValue toJSONValueImpl(F)(Complex!F dst)
    if(isFloatingPoint!F)
    out(result){
        assert(result.type == JSONType.object);
    }
    do{
        JSONValue json = JSONValue(string[string].init);
        json["re"] = dst.re;
        json["im"] = dst.im;
        return json;
    }


    ///
    JSONValue toJSONValueImpl(R)(R range)
    if(isInputRange!R && !isSomeString!R)
    out(result){
        assert(result.type == JSONType.array);
    }
    do{
        auto app = appender!(JSONValue[])();

        app.put(range.map!(a => toJSONValue(a)));

        return JSONValue(app.data);
    }


    ///
    JSONValue toJSONValueImpl(AA)(AA aa)
    if(isAssociativeArray!AA)
    out(result){
        assert(result.type == JSONType.object);
    }
    do{
        JSONValue[string] dst;
        foreach(k, v; aa){
            static if(is(typeof(k) : string))
                dst[k] = toJSONValue(v);
            else
                dst[k.to!string()] = toJSONValue(v);
        }

        return JSONValue(dst);
    }


    /// 
    JSONValue toJSONValueImpl(T)(T v)
    if(is(T == Variant))
    {
        if(!v.hasValue)
            return JSONValue(null);
        else
        {
            bool bSuccess;
            JSONValue jv;

            alias TT = TypeTuple!(Variant, Variant[], Variant[string], Variant[Variant],
                                  byte, ubyte, short, ushort,
                                  int, uint, long, ulong,
                                  float, double, real,
                                  string, wstring, dstring,
                                  bool, typeof(null));

            foreach(U; TT){
                if(auto p = v.peek!U){
                    jv = toJSONValue(*p);
                    bSuccess =  true;
                    break;
                }
            }

            if(!bSuccess)
                jv = JSONValue(v.to!string());

            return jv;
        }
    }


    private string createFromJSONValueExceptionMsg(T)(JSONValue json)
    {
        return "cannot convert to '" ~ T.stringof ~ "' from " ~ `"` ~ toJSON(json) ~ `"`;
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(is(T == JSONValue))
    {
        dst = json;
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(is(T == typeof(null)))
    {
        enforce!JSONException(json.type == JSONType.null_, createFromJSONValueExceptionMsg!T(json));
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isSomeString!T)
    {
        enforce!JSONException(json.type == JSONType.string, createFromJSONValueExceptionMsg!T(json));
        dst = json.str.to!T;
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isIntegral!T && isUnsigned!T)
    {
        enforce!JSONException(json.type == JSONType.uinteger || json.type == JSONType.integer, createFromJSONValueExceptionMsg!T(json));

        if(json.type == JSONType.uinteger)
            dst = json.uinteger.to!T();
        else
            dst = json.integer.to!T();
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isIntegral!T && isSigned!T)
    {
        enforce!JSONException(json.type == JSONType.integer || json.type == JSONType.uinteger, createFromJSONValueExceptionMsg!T(json));

        if(json.type == JSONType.integer)
            dst = json.integer.to!T();
        else
            dst = json.uinteger.to!T();
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(is(T == bool))
    {
        enforce!JSONException(json.type == JSONType.true_ || json.type == JSONType.false_, createFromJSONValueExceptionMsg!T(json));
        dst = json.type == JSONType.true_;
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isFloatingPoint!T)
    {
        enforce!JSONException(json.type == JSONType.float_
                             || json.type == JSONType.integer
                             || json.type == JSONType.uinteger, createFromJSONValueExceptionMsg!T(json));
        
        if(json.type == JSONType.float_)
            dst = json.floating;
        else if(json.type == JSONType.integer)
            dst = json.integer;
        else
            dst = json.uinteger;
    }


    ///
    void fromJSONValueImpl(F)(JSONValue json, ref Complex!F dst)
    if(isFloatingPoint!F)
    {
        enforce!JSONException(json.type == JSONType.object);
        enforce!JSONException("re" in json.object && "im" in json.object);
        enforce!JSONException(json["re"].type == JSONType.float_ && json["im"].type == JSONType.float_);

        dst.re = json["re"];
        dst.im = json["im"];
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isArray!T && !isSomeString!T)
    {
        enforce!JSONException(json.type == JSONType.array, createFromJSONValueExceptionMsg!T(json));

        T data = new T(json.array.length);

        foreach(i, e; json.array){
            typeof(data[i]) elem;
            fromJSONValue(e, elem);
            data[i] = elem;
        }

        dst = data;
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isInputRange!T && isOutputRange!(T, Unqual!(ElementType!T)) && !isArray!T)
    {
        enforce!JSONException(json.type == JSONType.array, createFromJSONValueExceptionMsg!T(json));

        foreach(e; json.array){
            alias Elem = Unqual!(ElementType!T);
            Elem elem;
            fromJSONValue(e, elem);
            dst.put(elem);
        }
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(isAssociativeArray!(T))
    {
        enforce!JSONException(json.type == JSONType.object, createFromJSONValueExceptionMsg!T(json));

        alias V = typeof(T.init.values[0]);
        alias K = typeof(T.init.keys[0]);

        foreach(k, v; json.object){
            V elem;
            fromJSONValue(v, elem);
            dst[k.to!K] = elem;
        }
    }


    ///
    void fromJSONValueImpl(T)(JSONValue json, ref T dst)
    if(is(T == Variant))
    {
        static void impl(T)(JSONValue json, ref Variant dst)
        {
            T v;
            fromJSONValue(json, v);
            dst = v;
        }


        final switch(json.type)
        {
          case JSONType.string:
            impl!string(json, dst);
            return;

          case JSONType.integer:
            impl!long(json, dst);
            return;

          case JSONType.uinteger:
            impl!ulong(json, dst);
            return;

          case JSONType.float_:
            impl!real(json, dst);
            return;

          case JSONType.object:
            impl!(Variant[string])(json, dst);
            return;

          case JSONType.array:
            impl!(Variant[])(json, dst);
            return;

          case JSONType.true_, JSONType.false_:
            impl!(bool)(json, dst);
            return;

          case JSONType.null_:
            impl!(typeof(null))(json, dst);
            return;
        }
    }
}


alias DefaultJSONEnv = JSONEnv!(null);


/**
任意のユーザー定義型をJSONでのObjectに変換する際に便利なテンプレート
*/
mixin template JSONObject(fields...)
if(fields.length && fields.length % 2 == 0)
{
    JSONValue toJSONValueImpl() @property
    {
        JSONValue[string] aa;

        foreach(i; _StaticIota!(0, fields.length))
        {
            static if(i % 2 == 0)
            {
                static assert(is(typeof(fields[i]) == string));

                static if(is(typeof(mixin(fields[i+1]))))
                    aa[fields[i]] = toJSONValue(mixin(fields[i+1]));
                else
                    aa[fields[i]] = toJSONValue(fields[i+1]);
            }
        }

        return JSONValue(aa);
    }


    private ref typeof(this) fromJSONValueImpl_(JSONValue jv) return
    {
        foreach(i; _StaticIota!(0, fields.length))
        {
            static if(i % 2 == 0)
            {
                static if(is(typeof(&(fields[i+1]())) == U*, U))
                    fromJSONValue(jv.object[fields[i]], fields[i+1]());
                else static if(is(typeof(&(fields[i+1])) == U*, U))
                    fromJSONValue(jv.object[fields[i]], fields[i+1]);
                else static if(is(typeof(&(mixin(fields[i+1]))) == U*, U) &&
                               !is(typeof(&(mixin(fields[i+1]))) == V function(W), V, W...))
                {
                    fromJSONValue(jv.object[fields[i]], mixin(fields[i+1]));
                }
                else
                {
                    static if(is(typeof(fields[i+1]())))    // property
                        alias X = typeof({auto x = fields[i+1](); return x;}());
                    else static if(is(typeof(mixin(fields[i+1])())))    // property(string)
                        alias X = typeof(mixin(fields[i+1])());
                    else
                        alias X = ParameterTypeTuple!(fields[i+1])[0];

                    X x;
                    fromJSONValue(jv.object[fields[i]], x);

                    static if(is(typeof(mixin(fields[i+1]))))
                        mixin(fields[i+1]) = x;
                    else
                        fields[i+1](x);
                }
            }
        }


        return this;
    }


    static typeof(this) fromJSONValueImpl(JSONValue jv)
    {
        typeof(this) dst;
        dst.fromJSONValueImpl_(jv);
        return dst;
    }
}


///
unittest{
    // Custom JSON Convertor
    // ユーザーは、任意の型をJSONへ変換するための変換器を定義できる
    static struct CustomJSONConvertor
    {
        static { mixin JSONEnv!null _dummy; }

        static JSONValue toJSONValueImpl(string str)
        {
            return _dummy.toJSONValueImpl("Custom Convertor-String : " ~ str);
        }


        static void fromJSONValueImpl(JSONValue json, ref string str)
        {
            assert(json.type == JSONType.string, "Error");
            str = json.str.find(" : ").drop(3);
        }
    }


    // グローバルなsetterとgetterだと仮定
    static struct Foo
    {
        int gloF() @property { return _gloF; }
        void gloF(int a) @property { _gloF = a; }

        static int _gloF = 12345;
    }


    // グローバル変数だと仮定
    static string gloG = "global variable";


    // JSONへ変換したり、JSONから変換する対象のオブジェクト
    static struct S1
    {
        // JSON Convertorの定義
        // 通常はモジュールで定義すればよい
        static {
            mixin JSONEnv!CustomJSONConvertor;
        }


        // refで返すプロパティ
        ref real flt() @property return
        {
            return _flt;
        }


        // getter
        int[] arr() @property
        {
            return _arr;
        }

        // setter
        void arr(int[] arr) @property
        {
            _arr = arr;
        }


        // JSONでのオブジェクトの定義
        mixin JSONObject!("intA", a,                // メンバ変数
                          "strB", "b",              // メンバ変数(文字列)
                          "fltC", flt,              // refを返すメンバ関数(プロパティ)
                          "arrD", "arr",            // setter, getter
                          "aasE", "aas",            // static setter, getter
                          "gloF", "Foo.init.gloF",  // global setter, getter
                          "gloG", gloG,             // グローバル変数など外部スコープの変数
                          );


      private:
        int a;
        string b;
        real _flt;
        int[] _arr;


      static:

        // staticなgetter
        int[int] aas() @property
        {
            return [1 : 2, 2 : 3, 3 : 4, 4 : 5];
        }


        // staticなsetter
        void aas(int[int] aa) @property
        {}
    }

    auto s1 = S1(12, "foo", 2.0, [1, 2, 3, 4]);
    auto jv = S1.toJSONValue(s1);

    auto jvtext = parseJSON(`{"gloF":12345,"strB":"Custom Convertor-String : foo","fltC":2.0,"gloG":"Custom Convertor-String : global variable","intA":12,"aasE":{"4":5,"1":2,"2":3,"3":4},"arrD":[1,2,3,4]}`);
    assert(toJSON(jv) == toJSON(jvtext));

    auto s2 = S1.fromJSONValueImpl(jv);
    auto s3 = S1.fromJSONValueImpl(jvtext);

    assert(s1 == s2);
    assert(s2 == s3);

    assert(Foo.init.gloF == 12345);
    assert(gloG == "global variable");
}


/**
JSON Objectを構築します
*/
struct JSONObjectType(alias jsonEnv, fields...)
if(fields.length && isValidJSONObjectTypeFields!fields)
{
    static {
        mixin JSONEnv!jsonEnv;
    }

    mixin(genStructField);
    mixin(`mixin JSONObject!(` ~ genMixinFields ~ `);`);

  private:
  static:
    string genStructField()
    {
        string dst;

        foreach(i; _StaticIota!(0, fields.length))
        {
            static if(i % 2 == 0)
            {
                dst ~= fields[i].stringof;  // type
                dst ~= " ";
                dst ~= fields[i+1];         // identifier
                dst ~= ";\n";
            }
        }

        return dst;
    }


    string genMixinFields()
    {
        string dst;

        foreach(i; _StaticIota!(0, fields.length))
        {
            static if(i % 2 == 0)
            {
                dst ~= `"` ~ fields[i+1] ~ `",`;     // tag
                dst ~= `"` ~ fields[i+1] ~ `",`;     // identifier
            }
        }

        return dst;
    }
}


unittest{
    JSONObjectType!(null,
                    int,        "intA",
                    string,     "strB",
                    float,      "fltC") jot;

    jot.intA = 12;
    jot.strB = "foo-bar";
    jot.fltC = 12.5;

    auto x = jot.toJSONValueImpl();
    auto y = parseJSON(`{"strB":"foo-bar","fltC":12.5,"intA":12}`);
    assert(toJSON(x) == toJSON(y));

    assert(typeof(jot).fromJSONValueImpl(x) == jot);
    assert(typeof(jot).fromJSONValueImpl(y) == jot);
}


private
template isValidJSONObjectTypeFields(fields...)
{
    //template isValidImpl(T, string name) { enum isValidImpl = true; }
    enum isValidImpl_(T, string name) = name.length;
    enum isValidImpl(X...) = is(typeof({static assert(isValidImpl_!X);}));

    static if(fields.length)
    {
        static if(fields.length % 2 == 0)
            enum isValidJSONObjectTypeFields = isValidImpl!(fields[0 .. 2])
                                            && isValidJSONObjectTypeFields!(fields[2 .. $]);
        else
            enum isValidJSONObjectTypeFields = false;
    }
    else
        enum isValidJSONObjectTypeFields = true;
}
