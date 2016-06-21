module dffdd.blockdiagram.filter;


import std.range;
import std.traits;


struct FIRFilter
{
    static
    auto makeBlock(R, C)(R r, C[] coefs)
    if(isInputRange!R)
    {
        return SimpleFIRFilter!(R, C)(r, coefs);
    }


    static
    auto makeBlock(R, C)(R r, C coef1)
    if(isInputRange!R && is(typeof((R r, C c){ auto e = r.front * c; })))
    {
        return Simple1TapFIRFilter!(R, C)(r, coef1);
    }



    static struct SimpleFIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);


        this(R r, C[] coefs)
        {
            _r = r;
            _coefs = coefs.dup;
            _regs = (new E[coefs.length]).cycle();
            /*
            foreach(i; 0 .. coefs.length-1){
                if(_r.empty){
                    _empty = true;
                    return;
                }

                _regs.popFront();
                _regs[_coefs.length - 1] = _r.front;
                _r.popFront();
            }*/
            foreach(i; 0 .. coefs.length)
                _regs[i] = cast(F)0;

            this.popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }

        void popFront()
        {
            _regs.popFront();

            if(_r.empty){
                _empty = true;
                return;
            }
            _regs[_coefs.length - 1] = _r.front;
            _r.popFront();
            _front = cast(F)0;
            foreach(i; 0 .. _coefs.length)
                _front += _coefs[$-1-i] * _regs[i];
        }


      private:
        R _r;
        const(C)[] _coefs;
        typeof(cycle(new E[1])) _regs;
        F _front;
        bool _empty;
    }


    static struct Simple1TapFIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);


        this(R)(R r, C coef1)
        {
            _r = r;
            _coef1 = coef1;
            _x0 = cast(E)0;
            _x1 = cast(E)0;
            popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }


        void popFront()
        {
            _x1 = _x0;

            if(_r.empty){
                _empty = true;
                return;
            }

            _x0 = _r.front;
            _r.popFront();

            _front = _x0 + _x1*_coef1;
        }


      private:
        R _r;
        C _coef1;
        E _x0, _x1;
        F _front;
        bool _empty;
    }
}


struct IIRFilter
{
    static
    auto makeBlock(R, C)(R r, C coef1)
    if(isInputRange!R)
    {
        return Simple1TapIIRFilter!(R, C)(r, coef1);
    }


    static struct Simple1TapIIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);

        this(R)(R r, C coef1)
        {
            _r = r;
            _coef1 = coef1;
            _front = cast(F)0;
            popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }


        void popFront()
        {
            if(_r.empty){
                _empty = true;
                return;
            }

            _front = _r.front + _coef1 * _front;
            _r.popFront();
        }

      private:
        R _r;
        C _coef1;
        F _front;
        bool _empty;
    }
}
