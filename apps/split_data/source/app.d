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

import carbon.stream;

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
    //auto rawStream = rawFileStream!ubyte(input);
    File ifile = File(input);
    ubyte[] buf = new ubyte[1024*1024];

    foreach(ubyte[] chunk; ifile.byChunk(buf)){
        if(offset){
            auto m = min(chunk.length, offset);
            offset -= m;
            chunk = chunk[m .. $];
        }

        if(total && chunk.length){
            auto m = min(total, chunk.length);
            total -= m;
            ofile.rawWrite(chunk[0 .. m]);
        }

        if(!total) break;
    }

    //while(offset) { offset -= rawStream.read(buf[0 .. min($, offset)]).length; writeln(offset); }
    //while(total){
    //    auto rs = rawStream.read(buf[0 .. min($, total)]);
    //    total -= rs.length;
    //    ofile.rawWrite(rs);
    //}
}
