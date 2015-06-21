import std.stdio;
import std.getopt;
import std.conv;
import std.format;
import std.string;
import std.algorithm;

import dffdd.gps.code;

enum GPSCode { L1CA, L2CM, L2CL, MIX }

void main(string[] args)
{
    GPSCode code = GPSCode.L1CA;
    uint prn = 1;
    string filename;

    auto opt = getopt(
        args,
        "prn", &prn,
        "code", &code,
        "output", &filename
    );

    if(opt.helpWanted)
    {
        defaultGetoptPrinter("Generate GPS Baseband code such as L1CA, L2CM and L2CL.",
            opt.options);
    }

    if(filename is null) filename = format("%s.dat", code);

    final switch(code)
    {
      case GPSCode.L1CA:
        auto sig = L1CACode.repeatStream(prn);
        outputToFile(sig, filename, L1CACode.codeLength);
        break;

      case GPSCode.L2CM:
        auto sig = L2CMCode.repeatStream(prn);
        outputToFile(sig, filename, L2CMCode.codeLength);
        break;

      case GPSCode.L2CL:
        auto sig = L2CLCode.repeatStream(prn);
        outputToFile(sig, filename, L2CLCode.codeLength);
        break;

      case GPSCode.MIX:
        outputMIXToFile(prn, filename, L2CLCode.codeLength*3);
    }
}


void outputToFile(S)(S signal, string filename, ptrdiff_t numSamples)
{
    auto file = File(filename, "w");
    auto bits = new byte[1024*8];

    while(numSamples > 0)
    {
        auto res = signal.readOp!"a=(b==1?0:1)"(bits[0 .. min(numSamples*8, $)]);
        assert(res.length % 8 == 0);
        immutable nbyte = res.length / 8;
        assert(nbyte <= numSamples);

        foreach(idx; 0 .. nbyte)
        {
            auto i = idx * 8;
            res[idx] = 0;
            foreach(j; 0 .. 8)
                res[idx] += res[i+j] << j;
        }

        file.rawWrite(res[0 .. nbyte]);
        numSamples -= nbyte*8;
    }
}


void outputMIXToFile(uint prn, string filename, ptrdiff_t numSamples)
{
    auto file = File(filename, "w");
    auto l1ca = L1CACode.repeatStream(prn);
    auto l2cm = L2CMCode.repeatStream(prn);
    auto l2cl = L2CLCode.repeatStream(prn);

    while(numSamples > 0)
    {
        ubyte[8][3] bits;
        assert(l1ca.readOp!"a=(b==1?0:1)"(bits[0][]).length == 8);
        assert(l2cm.readOp!"a=(b==1?0:1)"(bits[1][]).length == 8);
        assert(l2cl.readOp!"a=(b==1?0:1)"(bits[2][]).length == 8);

        uint uint3;
        foreach(i; 0 .. 3)
            foreach(j; 0 .. 8)
                uint3 += bits[i][j] << (i + j*3);

        ubyte[] data = (cast(ubyte*)&uint3)[0 .. 3];
        file.rawWrite(data);
        numSamples -= 3*8;
    }
}
