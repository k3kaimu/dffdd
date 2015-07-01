import std.stdio;
import std.getopt;
import std.range;
import std.algorithm;
import std.array;
import std.format;
import std.string;
import std.numeric;
import std.math;
import std.complex;

void main(string[] args)
{
    size_t offset, total;
    string input = "input.dat", output = "output.dat";

    getopt(args,
        "offset", &offset,
        "total", &total,
        "input", &input,
        "output", &output);

    File ofile = File(output, "w");
    File ifile = File(input);
    cfloat[] buf = new cfloat[1024*1024];

    ifile.seek(offset * cfloat.sizeof);
    while(total){
        auto rs = ifile.rawRead(buf[0 .. min($, total)]);
        total -= rs.length;

        foreach(e; rs)
            ofile.writefln("%s,%s,", e.re, e.im);
    }
}
