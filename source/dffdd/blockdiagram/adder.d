module dffdd.blockdiagram.adder;

import std.algorithm;
import std.range;


auto add(R1, R2)(R1 r1, R2 r2)
{
    return r1.zip(r2).map!"a[0]+a[1]";
}
