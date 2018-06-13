module dffdd.utils.optimize;

__EOF__

import std.math;
import std.complex;
import std.typecons;
import libnlopt;

/**
目的関数f(x)を最小にするxを傾きdf(x)/dxを用いて求めます．
*/
C newtonMethod(Fn, Gn, C)(Fn fn, Gn dfn, C initX, real sqtol = 1E-4)
{
    C x = initX;
    real err = 1;
    do{
        C fx = fn(x);
        C dfx = dfn(x);
        x -= fx / dfx;

        static if(is(typeof(x.re)))
            err = fx.sqAbs;
        else
            err = fx^^2;
    }while(err > sqtol);

    return x;
}


unittest
{
    import std.stdio;
    import std.complex;
    alias C = Complex!real;

    // x - x|x|^2 = -1-i となるようなxを求めたい
    C ans = newtonMethod(
        (C x) => x - x * x.sqAbs - C(-1, -1),
        (C x) => 1 - 2 * x.sqAbs,
        C(0, 0)
    );

    assert(sqAbs(ans - ans * ans.sqAbs - C(-1, -1)) < 1E-4);
}


/**
*/
struct NLopt
{
    this(nlopt_algorithm alg, uint dim)
    {
        _nlopt = RefCounted!NLoptInstance(alg, dim);
    }


    void setLowerBounds(in double[] lb)
    {
        nlopt_set_lower_bounds(_nlopt.handle, lb.ptr);
    }


    void setLowerBounds1(double lb)
    {
        nlopt_set_lower_bounds(_nlopt.handle, lb);
    }


    void setUpperBounds(in double[] ub)
    {
        nlopt_set_upper_bounds(_nlopt.handle, ub.ptr);
    }


    void setUpperBounds1(double ub)
    {
        nlopt_set_upper_bounds(_nlopt.handle, ub);
    }


    void removeInequalityConstraints()
    {
        nlopt_remove_inequality_constraints(_nlopt.handle);
        _nlopt.allocedIneqConList.length = 0;
    }


    void addInequalityConstraint(Fn)(Fn fn, double tol)
    {
        Fn* fptr = new Fn(fn);
        allocedIneqConList ~= cast(void*)fptr;
        nlopt_add_inequality_constraint(_nlopt.handle, &singlefunction!Fn, cast(void*)fptr, tol);
    }


    void addInequalityMConstraint(Fn)(uint m, Fn fn, in double[] tols)
    {
        Fn* fptr = new Fn(fn);
        allocedIneqConList ~= cast(void*)fptr;
        nlopt_add_inequality_mconstraint(_nlopt.handle, m, &multifunction!Fn, cast(void*)fptr, tol);
    }


    void removeEqualityConstraits()
    {
        nlopt_remove_equality_constraints(_nlopt.handle);
        allocedEqConList.length = 0;
    }


    void addEqualityConstraint(Fn)(Fn fn, double tol)
    {
        Fn* fptr = new Fn(fn);
        allocedEqConList ~= cast(void*)fptr;
        nlopt_add_equality_constraint(_nlopt.handle, &singlefunction!Fn, cast(void*)fptr, tol);
    }


    void addEqualityMConstraint(Fn)(Fn fn, in double[] tols)
    {
        Fn* fptr = new Fn(fn);
        allocedEqConList ~= cast(void*)fptr;
        nlopt_add_equality_mconstraint(_nlopt.handle, &multifunction!Fn, cast(void*)fptr, tol);
    }


  private:
    RefCounted!NLoptInstance _nlopt;


    static struct NLoptInstance
    {
        this(nlopt_algorithm alg, uint dim)
        {
            handle = nlopt_create(alg, dim);
        }


        ~this()
        {
            nlopt_destroy(handle);
        }


        nlopt_opt handle;
        void*[] allocedIneqConList;
        void*[] allocedEqConList;
    }


    static
    double singlefunction(F)(uint dim, const double* x, double* grad, void* data)
    {
        auto func = *cast(F*)data;
    }
}
