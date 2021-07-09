module dffdd.gps.code;

import std.algorithm;
import std.array;
import std.range;


import carbon.stream;

import dffdd.utils.binary;


/**
Gold系列のコードを生成します
*/
struct L1CACode
{
    enum size_t codeLength = 1023;

    shared static this()
    {
        byte[10] r1, r2;
        r1[] = cast(byte)-1;
        r2[] = cast(byte)-1;

        int[codeLength] g1, g2;
        foreach(i; 0 .. codeLength){
            g1[i] = r1[9];
            g2[i] = r2[9];

            byte c1 = cast(byte)(r1[2] * r1[9]),
                 c2 = cast(byte)(r2[1] * r2[2] * r2[5] * r2[7] * r2[8] * r2[9]);

            foreach_reverse(j; 1 .. 10){
                r1[j] = r1[j-1];
                r2[j] = r2[j-1];
            }

            r1[0] = c1;
            r2[0] = c2;
        }

        _g1[] = g1[];
        _g2[] = g2[];
    }


  private static
  {
    shared immutable int[codeLength] _g1, _g2;

    shared immutable short[] _delay = cast(short[])[ /* G2 delay (chips) */
          5,   6,   7,   8,  17,  18, 139, 140, 141, 251,   /*   1- 10 */
        252, 254, 255, 256, 257, 258, 469, 470, 471, 472,   /*  11- 20 */
        473, 474, 509, 512, 513, 514, 515, 516, 859, 860,   /*  21- 30 */
        861, 862, 863, 950, 947, 948, 950,  67, 103,  91,   /*  31- 40 */
         19, 679, 225, 625, 946, 638, 161,1001, 554, 280,   /*  41- 50 */
        710, 709, 775, 864, 558, 220, 397,  55, 898, 759,   /*  51- 60 */
        367, 299,1018, 729, 695, 780, 801, 788, 732,  34,   /*  61- 70 */
        320, 327, 389, 407, 525, 405, 221, 761, 260, 326,   /*  71- 80 */
        955, 653, 699, 422, 188, 438, 959, 539, 879, 677,   /*  81- 90 */
        586, 153, 792, 814, 446, 264,1015, 278, 536, 819,   /*  91-100 */
        156, 957, 159, 712, 885, 461, 248, 713, 126, 807,   /* 101-110 */
        279, 122, 197, 693, 632, 771, 467, 647, 203, 145,   /* 111-120 */
        175,  52,  21, 237, 235, 886, 657, 634, 762, 355,   /* 121-130 */
       1012, 176, 603, 130, 359, 595,  68, 386, 797, 456,   /* 131-140 */
        499, 883, 307, 127, 211, 121, 118, 163, 628, 853,   /* 141-150 */
        484, 289, 811, 202,1021, 463, 568, 904, 670, 230,   /* 151-160 */
        911, 684, 309, 644, 932,  12, 314, 891, 212, 185,   /* 161-170 */
        675, 503, 150, 395, 345, 846, 798, 992, 357, 995,   /* 171-180 */
        877, 112, 144, 476, 193, 109, 445, 291,  87, 399,   /* 181-190 */
        292, 901, 339, 208, 711, 189, 263, 537, 663, 942,   /* 191-200 */
        173, 900,  30, 500, 935, 556, 373,  85, 652, 310    /* 201-210 */
    ];
  }


    static
    auto lutNCO(uint prn, real freq, real deltaT, real theta = 0)
    {
        byte[] bs = new byte[codeLength];
        auto l1ca = L1CACode(prn);
        l1ca.read(bs);

        return .lutNCO(bs, freq, deltaT, theta);
    }


    static
    auto repeatStream(uint prn)
    {
        byte[] bs = new byte[codeLength];
        auto l1ca = L1CACode(prn);
        l1ca.read(bs);

        return .repeatStream(bs);
    }


    this(uint prn)
    in{
        assert(prn > 0);
    }
    body{
        _prn = prn;
        _i = 0;
        _j = codeLength - _delay[prn-1];
    }

    byte front() const @property { return cast(byte)(-_g1[_i] * _g2[_j % $]); }
    void popFront(){ ++_i; ++_j; _i %= codeLength; _j %= codeLength; }
    enum bool empty = false;
    L1CACode save() const { return this; }
    byte opIndex(size_t i) const { return cast(byte)(-_g1[(_i + i)%$] * _g2[(_j + i)%$]); }
    struct OpDollar{} enum opDollar = OpDollar();
    L1CACode opSlice() const { return this; }
    L1CACode opSlice(size_t i, OpDollar) const { auto dst = this.save; dst._i = (_i + i) % codeLength; dst._j = (_j + i) % codeLength; return dst; }
    auto opSlice(size_t i, size_t j) const { return take(this[i .. $], j-i); }

    uint prn() const @property { return _prn; }
    void prn(uint p) @property { _prn = p; _j = (codeLength - _delay[p-1] + _i) % codeLength; }

    bool fetch() { return false; }

    E[] readOp(alias op, E)(E[] buf)
    {
        immutable len = buf.length;
        {
            auto p = buf.ptr;
            const end = () @trusted { return p + buf.length; }();
            size_t i = _i % codeLength,
                   j = _j % codeLength,
                   k = 0;

            while(p != end && k != codeLength){
                binaryFunExt!op(*p, cast(byte)(-_g1[i%$] * _g2[j%$]));
                ++i; ++j; ++k; ++p;
            }
        }

        if(len > codeLength){
            foreach(k; 1 .. len / codeLength){
                auto p = () @trusted { return buf.ptr + k*codeLength; }();
                const end = () @trusted { return p + codeLength; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }

            immutable rem = len % codeLength;
            buf[$ - rem .. $] = buf[0 .. rem];
            {
                auto p = () @trusted { return buf.ptr + len - rem; }();
                const end = () @trusted { return buf.ptr + len; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }
        }

        _i = (_i + buf.length) % codeLength;
        _j = (_j + buf.length) % codeLength;
        return buf;
    }


    T[] read(T)(T[] buf) @property { return readOp!""(buf); }


  private:
    uint _prn;
    size_t _i, _j;
}

unittest
{
    import std.stdio : writeln;

    auto l1ca1 = L1CACode(1);
    assert(equal(l1ca1.take(8), [1, 1, -1, -1, 1, -1, -1, -1]));

    auto l1ca2 = L1CACode(2);
    assert(equal(l1ca2.take(8), [1, 1, 1, -1, -1, 1, -1, -1]));

    // みちびき
    auto l1ca193 = L1CACode(193);
    assert(equal(l1ca193.take(8), [-1, 1, 1, 1, -1, 1, -1, 1]));

    immutable len = l1ca193.codeLength;

    auto buf = new byte[len * 3 + 10];
    assert(l1ca193.read(buf) == buf);

    assert(buf[0 .. len] == buf[len .. len*2]);
    assert(buf[len .. len*2] == buf[len*2 .. len*3]);
    assert(buf[0 .. 10] == buf[len*3 .. len*3+10]);
}



struct L2CMCode
{
    enum size_t codeLength = 10230;

    static immutable string[210] reg = [
      "742417664","756014035","002747144","066265724","601403471","703232733","124510070","617316361","047541621","733031046",   /*   1- 10 */
      "713512145","024437606","021264003","230655351","001314400","222021506","540264026","205521705","064022144","120161274",   /*  11- 20 */
      "044023533","724744327","045743577","741201660","700274134","010247261","713433445","737324162","311627434","710452007",   /*  21- 30 */
      "722462133","050172213","500653703","755077436","136717361","756675453","435506112","771353753","226107701","022025110",   /*  31- 40 */
      "402466344","752566114","702011164","041216771","047457275","266333164","713167356","060546335","355173035","617201036",   /*  41- 50 */
      "157465571","767360553","023127030","431343777","747317317","045706125","002744276","060036467","217744147","603340174",   /*  51- 60 */
      "326616775","063240065","111460621","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  61- 70 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  71- 80 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  81- 90 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  91-100 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 101-110 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 111-120 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 121-130 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 131-140 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 141-150 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","604055104","157065232",   /* 151-160 */
      "013305707","603552017","230461355","603653437","652346475","743107103","401521277","167335110","014013575","362051132",   /* 161-170 */
      "617753265","216363634","755561123","365304033","625025543","054420334","415473671","662364360","373446602","417564100",   /* 171-180 */
      "000526452","226631300","113752074","706134401","041352546","664630154","276524255","714720530","714051771","044526647",   /* 181-190 */
      "207164322","262120161","204244652","202133131","714351204","657127260","130567507","670517677","607275514","045413633",   /* 191-200 */
      "212645405","613700455","706202440","705056276","020373522","746013617","132720621","434015513","566721727","140633660"    /* 201-210 */
    ];


    static
    auto lutNCO(uint prn, real freq, real deltaT, real theta = 0)
    {
        auto cs = L2CMCode(prn);
        return .lutNCO(cs._code, freq, deltaT, theta);
    }


    static
    auto repeatStream(uint prn)
    {
        auto cs = L2CMCode(prn);
        return .repeatStream(cs._code);
    }


    this(uint prn)
    {
        import std.exception : assumeUnique;

        auto code = new int[codeLength];

        int[27] r, t;
        oct2bin(reg[prn-1], r[], 0, 0);     // initial state
        oct2bin("112225170", t[], 0, 0);    // polynomial coefficient

        foreach(i; 0 .. codeLength){
            int c = r[26];
            code[i] = -c;

            foreach(j; 0 .. 27) if(t[j] ==  -1) r[j] *= c;
            foreach_reverse(j; 1 .. 27) r[j] = r[j-1];
            r[0] = c;
        }

        _code = code.map!"cast(byte)a".array().idup;
        _prn = prn;
        _i = 0;
    }


    byte front() const @property { return cast(byte)(_code[_i]); }
    void popFront(){ ++_i; _i %= codeLength; }
    enum bool empty = false;
    L2CMCode save() const { return this; }
    byte opIndex(size_t i) const { return cast(byte)(_code[(_i + i)%$]); }
    struct OpDollar{} enum opDollar = OpDollar();
    L2CMCode opSlice() const { return this; }
    L2CMCode opSlice(size_t i, OpDollar) const { auto dst = this.save; dst._i = (_i + i) % codeLength; return dst; }
    auto opSlice(size_t i, size_t j) const { return take(this[i .. $], j-i); }

    uint prn() const @property { return _prn; }

    bool fetch() { return false; }

    E[] readOp(alias op, E)(E[] buf)
    {
        immutable len = buf.length;
        {
            auto p = buf.ptr;
            const end = () @trusted { return p + buf.length; }();
            size_t i = _i % codeLength,
                   k = 0;

            while(p != end && k != codeLength){
                binaryFunExt!op(*p, _code[i]);
                ++i; ++k; ++p;
            }
        }

        if(len > codeLength){
            foreach(k; 1 .. len / codeLength){
                auto p = () @trusted { return buf.ptr + k*codeLength; }();
                const end = () @trusted { return p + codeLength; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }

            immutable rem = len % codeLength;
            buf[$ - rem .. $] = buf[0 .. rem];
            {
                auto p = () @trusted { return buf.ptr + len - rem; }();
                const end = () @trusted { return buf.ptr + len; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }
        }

        _i = (_i + buf.length) % codeLength;
        return buf;
    }


    T[] read(T)(T[] buf) @property { return readOp!""(buf); }


  private:
    uint _prn;
    size_t _i;
    immutable(byte)[] _code;
}

unittest
{
    import std.stdio : writeln;

    auto l2cm1 = L2CMCode(1);
    assert(equal(l2cm1.take(8), [-1, -1, 1, -1, 1, -1, 1, 1]));

    auto l2cm2 = L2CMCode(2);
    assert(equal(l2cm2.take(8), [1, -1, 1, -1, -1, -1, -1, 1]));

    // みちびき
    auto l2cm193 = L2CMCode(193);
    assert(equal(l2cm193.take(8), [-1, 1, -1, 1, 1, -1, -1, -1]));

    immutable len = l2cm193.codeLength;

    auto buf = new byte[len * 3 + 10];
    assert(l2cm193.read(buf) == buf);

    assert(buf[0 .. len] == buf[len .. len*2]);
    assert(buf[len .. len*2] == buf[len*2 .. len*3]);
    assert(buf[0 .. 10] == buf[len*3 .. len*3+10]);
}



struct L2CLCode
{
    enum codeLength = 767250;

    static immutable string[] reg = [ /* Initial Register (Octal) */
      "624145772","506610362","220360016","710406104","001143345","053023326","652521276","206124777","015563374","561522076",   /*   1- 10 */
      "023163525","117776450","606516355","003037343","046515565","671511621","605402220","002576207","525163451","266527765",   /*  11- 20 */
      "006760703","501474556","743747443","615534726","763621420","720727474","700521043","222567263","132765304","746332245",   /*  21- 30 */
      "102300466","255231716","437661701","717047302","222614207","561123307","240713073","101232630","132525726","315216367",   /*  31- 40 */
      "377046065","655351360","435776513","744242321","024346717","562646415","731455342","723352536","000013134","011566642",   /*  41- 50 */
      "475432222","463506741","617127534","026050332","733774235","751477772","417631550","052247456","560404163","417751005",   /*  51- 60 */
      "004302173","715005045","001154457","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  61- 70 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  71- 80 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  81- 90 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /*  91-100 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 101-110 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 111-120 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 121-130 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 131-140 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000",   /* 141-150 */
      "000000000","000000000","000000000","000000000","000000000","000000000","000000000","000000000","605253024","063314262",   /* 151-160 */
      "066073422","737276117","737243704","067557532","227354537","704765502","044746712","720535263","733541364","270060042",   /* 161-170 */
      "737176640","133776704","005645427","704321074","137740372","056375464","704374004","216320123","011322115","761050112",   /* 171-180 */
      "725304036","721320336","443462103","510466244","745522652","373417061","225526762","047614504","034730440","453073141",   /* 181-190 */
      "533654510","377016461","235525312","507056307","221720061","520470122","603764120","145604016","051237167","033326347",   /* 191-200 */
      "534627074","645230164","000171400","022715417","135471311","137422057","714426456","640724672","501254540","513322453"    /* 201-210 */
    ];


    static
    auto lutNCO(uint prn, real freq, real deltaT, real theta = 0)
    {
        auto cs = L2CLCode(prn);
        return .lutNCO(cs._code, freq, deltaT, theta);
    }


    static
    auto repeatStream(uint prn)
    {
        auto cs = L2CLCode(prn);
        return .repeatStream(cs._code);
    }


    this(uint prn)
    {
        import std.exception : assumeUnique;

        auto code = new byte[codeLength];

        byte[27] r, t;
        oct2bin(reg[prn-1], r[], 0, 0);     // initial state
        oct2bin("112225170", t[], 0, 0);    // polynomial coefficient

        foreach(i; 0 .. codeLength){
            byte c = r[26];
            code[i] = cast(byte)(-cast(int)c);

            foreach(j; 0 .. 27) if(t[j] ==  -1) r[j] *= c;
            foreach_reverse(j; 1 .. 27) r[j] = r[j-1];
            r[0] = c;
        }

        _code = assumeUnique(code);
        _prn = prn;
        _i = 0;
    }


    byte front() const @property { return _code[_i]; }
    void popFront(){ ++_i; _i %= codeLength; }
    enum bool empty = false;
    L2CLCode save() const { return this; }
    byte opIndex(size_t i) const { return _code[(_i + i)%$]; }
    struct OpDollar{} enum opDollar = OpDollar();
    L2CLCode opSlice() const { return this; }
    L2CLCode opSlice(size_t i, OpDollar) const { auto dst = this.save; dst._i = (_i + i) % codeLength; return dst; }
    auto opSlice(size_t i, size_t j) const { return take(this[i .. $], j-i); }

    uint prn() const @property { return _prn; }

    bool fetch() { return false; }

    E[] readOp(alias op, E)(E[] buf)
    {
        immutable len = buf.length;
        {
            auto p = buf.ptr;
            const end = () @trusted { return p + buf.length; }();
            size_t i = _i % codeLength,
                   k = 0;

            while(p != end && k != codeLength){
                binaryFunExt!op(*p, _code[i]);
                ++i; ++k; ++p;
            }
        }

        if(len > codeLength){
            foreach(k; 1 .. len / codeLength){
                auto p = () @trusted { return buf.ptr + k*codeLength; }();
                const end = () @trusted { return p + codeLength; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }

            immutable rem = len % codeLength;
            buf[$ - rem .. $] = buf[0 .. rem];
            {
                auto p = () @trusted { return buf.ptr + len - rem; }();
                const end = () @trusted { return buf.ptr + len; }();
                const(E)* q = buf.ptr;

                while(p != end) { binaryFunExt!op(*p, *q); ++p; ++q; }
            }
        }

        _i = (_i + buf.length) % codeLength;
        return buf;
    }


    T[] read(T)(T[] buf) @property { return readOp!""(buf); }


  private:
    uint _prn;
    size_t _i;
    immutable(byte)[] _code;
}

unittest
{
    import std.stdio : writeln;

    auto l2cl1 = L2CLCode(1);
    assert(equal(l2cl1.take(8), [-1, 1, -1, 1, -1, -1, 1, 1]));

    auto l2cl2 = L2CLCode(2);
    assert(equal(l2cl2.take(8), [-1, 1, -1, -1, -1, -1, -1, -1]));

    // みちびき
    auto l2cl193 = L2CLCode(193);
    assert(equal(l2cl193.take(8), [-1, 1, -1, 1, 1, 1, 1, -1]));

    immutable len = l2cl193.codeLength;

    auto buf = new byte[len * 3 + 10];
    assert(l2cl193.read(buf) == buf);

    assert(buf[0 .. len] == buf[len .. len*2]);
    assert(buf[len .. len*2] == buf[len*2 .. len*3]);
    assert(buf[0 .. 10] == buf[len*3 .. len*3+10]);
}
