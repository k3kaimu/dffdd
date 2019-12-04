module dffdd.utils.range;

import std.range;
import std.traits;


E[] consumeFill(R, E)(ref R range, E[] buffer)
if(isAssignable!(E, ElementType!R))
{
    size_t count = 0;
    foreach(ref e; buffer) {
        if(range.empty) break;

        e = range.front;
        range.popFront();
        ++count;
    }


    return buffer[0 .. count];
}
