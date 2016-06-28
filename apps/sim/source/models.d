module models;


template Model(alias Constant)
{
    alias BasisFunctions = AliasSeq!(x => x,
                                  x => x.conj,
                                  x => x * (x.re^^2 + x.im^^2),
                                  x => x.conj * (x.re^^2 + x.im^^2),
                                  x => x * (x.re^^2 + x.im^^2)^^2,
                                  x => x.conj * (x.re^^2 + x.im^^2)^^2,
                                  x => x * (x.re^^2 + x.im^^2)^^4,
                                  x => x.conj * (x.re^^2 + x.im^^2)^^4,
                                  );


    auto randomBits(uint prn) {
        import dffdd.gps.code;
        return L2CLCode(prn).map!"cast(ubyte)(a < 0 ? 1 : 0)";
    }


    auto siRandomBits() { return randomBits(1); }
    auto desiredRandomBits() { return randomBits(193); }


    auto modOFDM(uint nOS)
    {
        return chainedMod(
            dffdd.mod.qam.QAM(Constant.QAM.arity),
            dffdd.mod.ofdm.OFDM(Constant.OFDM.numOfFFT, Constant.OFDM.numOfCP, Constant.OFDM.numOfSubcarrier, nOS),
        );
    }


    auto connectToModulator(R, Mod)(R r, Mod modObj)
    {
        return r
        .splitN(modObj.symInputLength)
        .tmap!(reverseArgs!mod, [0])(modObj)
        .joiner()
        .toWrappedRange
        ;
    }


    auto connectToDemodulator(R, Mod)(R r, Mod modObj/*, Bits bits*/)
    {
        return r
        .map!"cast(cfloat)a"
        .splitN(modObj.symOutputLength)
        //.connectToOFDMEqualizer(modObj, bits)
        .tmap!(reverseArgs!demod, [0])(modObj)
        .joiner()
        .toWrappedRange
        ;
    }


    auto connectToUpSampler(R)(R r)
    {
        return r
        .tmap!"a.repeat(b).array()"(Constant.UpDownSampler.scaleOfUpSampling).joiner()
        .connectTo!FIRFilter(Constant.UpDownSampler.decimationFIRFilterTaps)
        .toWrappedRange
        ;
    }



    auto connectToDownSampler(R)(R r)
    {
        return r
        .drop(Constant.UpDownSampler.numOfShift)
        .connectTo!FIRFilter(Constant.UpDownSampler.decimationFIRFilterTaps)
        .connectTo!(SimpleDecimator!(Constant.UpDownSampler.scaleOfUpSampling))
        .toWrappedRange
        ;
    }


    auto thermalNoise(real freq = Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling)
    {
        return ThermalNoise(freq, Constant.ThermalNoise.temperature);
    }


    auto connectToAWGN(R)(R r)
    {
        return r.add(thermalNoise()).toWrappedRange;
    }


    //auto connectToFlatRayleighFadingChannel(R)(R r)
    //{
    //    BoxMuller!Random bm;
    //    bm.seed(unpredictableSeed());
    //    bm.popFront();

    //    return r
    //    .zip(
    //        bm
    //        .tmap!"(a/SQRT2).repeat(b)"((Constant.OFDM.numOfFFT + Constant.OFDM.numOfCP)*Constant.OFDM.scaleOfUpSampling*Constant.UpDownSampler.scaleOfUpSampling*100)
    //        .joiner
    //    )
    //    .map!"a[0]*a[1]";
    //}


    //auto connectToOFDMEqualizer(R, Mod, Bits)(R r, Mod modObj, Bits bits)
    //{
    //    return OFDMEqualizer.makeBlock(
    //        r,
    //        modObj,
    //        Constant.OFDM.numOfFFT,
    //        Constant.OFDM.numOfCP,
    //        Constant.OFDM.numOfSubcarrier,
    //        Constant.OFDM.scaleOfUpSampling,
    //        bits).toWrappedRange;
    //}


    auto connectToTXIQMixer(R)(R r)
    {
        return r.connectTo!IQImbalance(0.dB, Constant.TXIQMixer.IIR.dB).toWrappedRange;
    }


    auto connectToRXIQMixer(R)(R r)
    {
        return r.connectTo!IQImbalance(0.dB, Constant.RXIQMixer.IIR.dB).toWrappedRange;
    }


    auto connectToPowerAmplifier(R)(R r)
    {
        auto v = (Constant.PA.TX_POWER - Constant.PA.GAIN).dBm;

        return r
        .connectTo!PowerControlAmplifier(v)
        .connectTo!PowerAmplifier(Constant.PA.GAIN.dB, Constant.PA.IIP3.dBm)
        .toWrappedRange;
    }


    auto connectToLNA(R)(R r)
    {
        return r
        .add(thermalNoise().tmap!"a*b"(Constant.LNA.NF.dB.gain - 1))
        .tmap!"a*b"(Constant.LNA.GAIN.dB.gain)
        .toWrappedRange;
    }


    auto connectToQuantizer(R)(R r)
    {
        return r
        .connectTo!PowerControlAmplifier((30 - Constant.OFDM.PAPR).dBm)
        .connectTo!SimpleQuantizer(Constant.Quantizer.numOfBits)
        .toWrappedRange;
    }


    auto txChain(R)(R r)
    {
        return r.connectToUpSampler().connectToTXIQMixer().connectToPowerAmplifier();
    }


    auto rxChain(R)(R r)
    {
        return r.connectToLNA().connectToRXIQMixer().connectToDownSampler().connectToQuantizer();
    }


    auto makeParallelHammersteinFilter(Mod)(Mod mod)
    {
        import dffdd.filter.diagonal;
        import dffdd.filter.lms;
        import dffdd.filter.ls;
        import dffdd.filter.rls;
        import dffdd.filter.mempoly;
        import dffdd.filter.polynomial;
        import dffdd.filter.orthogonalize;

        auto state = new MemoryPolynomialState!(Complex!float, 8, 4, 0, 0, false, true)(1);

        //writeln("Use");
        //auto adapter = new LMSAdapter!(typeof(state))(state, 0.001, 1024, 0.5);
        //auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
        auto adapter = lsAdapter(state, 10000);

        return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    }


    auto makeCascadeHammersteinFilter(Mod)(Mod mod)
    {
        import dffdd.filter.diagonal;
        import dffdd.filter.lms;
        import dffdd.filter.ls;
        import dffdd.filter.rls;
        import dffdd.filter.mempoly;
        import dffdd.filter.polynomial;
        import dffdd.filter.orthogonalize;
        import std.meta;
        import std.stdio;

        writeln("orthogonalizer start");

        auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BasisFunctions)();

        writeln("orthogonalizer setting up now...");

        {
            orthogonalizer.start();
            scope(exit)
                orthogonalizer.finish();

            writeln("orthogonalizer setting up now...");

            auto signal = randomBits(1).connectToModulator(mod);

            writeln("orthogonalizer setting up now...");
            //.put(orthogonalizer, signal.take(1024*400));
            Complex!float[] buf = new Complex!float[1024];
            foreach(i; 0 .. 1024){
                foreach(j; 0 .. 1024){
                    auto f = signal.front;
                    buf[j] = complex(f.re, f.im);
                    signal.popFront();
                }

                .put(orthogonalizer, buf);
            }

            writeln("orthogonalizer setting up now...");
        }

        Complex!float[][BasisFunctions.length] coefs;
        foreach(i, ref e; coefs){
            e = new Complex!float[BasisFunctions.length];
            orthogonalizer.getCoefs(i, e);
        }

        writeln("orthogonalizer is setup-ed");

        auto st1 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[0].dup);
        auto st12 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[0].dup);
        auto st1c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[1].dup);
        auto st12c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[1].dup);
        //auto st2 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        auto st3 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        auto st32 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        //auto st32 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        auto st3c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[3].dup);
        auto st32c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[3].dup);
        auto st5 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[4].dup);
        auto st52 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[4].dup);
        auto st5c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[5].dup);
        auto st52c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[5].dup);
        //auto st7 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[6].dup);
        //auto st7c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[7].dup);
        

        //auto st5 = (new FIRFilter!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[4].dup);
        //auto st7 = (new FIRFilter!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[5].dup);

        writeln("return filter");

        return serialFilter(
                //st12,
                //makeRLSAdapter(st12, 1 - 1E-4, 1E-7),
                //st12c,
                //makeRLSAdapter(st12c, 1 - 1E-4, 1E-7),
                //st32,
                //makeRLSAdapter(st32, 1 - 1E-4, 1E-7),
                //st32c,
                //makeRLSAdapter(st32c, 1 - 1E-4, 1E-7),
                //st52,
                //makeRLSAdapter(st52, 1 - 1E-4, 1E-7),
                //st52c,
                //makeRLSAdapter(st52c, 1 - 1E-4, 1E-7),
                //st12,
                //makeRLSAdapter(st12, 0.98, 1E-7),
                //st32,
                //makeRLSAdapter(st32, 0.98, 1E-7),
                //st1c_2,
                //makeRLSAdapter(st1c_2, 0.98, 1E-7),
                //st3_2,
                //makeRLSAdapter(st3_2, 0.92, 1E-7),
                //st3c_2,
                //makeRLSAdapter(st3c_2, 0.999, 1E-7),
                st1,
                lsAdapter(st1, 10000),
                //lmsAdapter(st1, 0.001, 1024, 0.5),
                //makeRLSAdapter(st1, 1 - 1E-4, 1E-7),
                st1c,
                lsAdapter(st1c, 10000),
                //lmsAdapter(st1c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st1c, 1 - 1E-4, 1E-7),
                //st2,
                //lsAdapter(st2, 10000),
                //lmsAdapter(st2, 0.002, 1024, 0.5),
                //makeRLSAdapter(st2, /*0.9997*/1, 1E-7),
                st3,
                lsAdapter(st3, 10000),
                //lmsAdapter(st3, 0.001, 1024, 0.5),
                //makeRLSAdapter(st3, 1 - 1E-4, 1E-7),
                st3c,
                lsAdapter(st3c, 10000),
                //lmsAdapter(st3c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st3c, 1 - 1E-4, 1E-7),
                //st4,
                //lsAdapter(st4, 10000),
                //lmsAdapter(st4, 0.002, 1024, 0.5),
                //makeRLSAdapter(st4, /*0.9997*/1, 1E-7),
                //st4a,
                //lmsAdapter(st4a, 0.0005, 1024, 0.5),
                st5,
                lsAdapter(st5, 10000),
                //lmsAdapter(st5, 0.001, 1024, 0.5),
                //makeRLSAdapter(st5, 1 - 1E-4, 1E-7),
                st5c,
                lsAdapter(st5c, 10000),
                //lmsAdapter(st5c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st5c, 1 - 1E-4, 1E-7),
                //st7,
                //lsAdapter(st7, 10000),
                //*********lmsAdapter(st7, 0.001, 1024, 0.5),
                //makeRLSAdapter(st7, 1 - 1E-6, 1E-7),
                //st12,
                //lsAdapter(st1, 10000),
                //lmsAdapter(st12, 0.010, 1024, 0.5),
                //makeRLSAdapter(st7, 0.9997, 1E-7),
                //*********st7c,
                //lsAdapter(st7c, 1000),
                //*********lmsAdapter(st7c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st7c, 1 - 1E-6, 1E-7),
                //st5,
                //lmsAdapter(st5, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st5, 1, 1E-7),
                //st3,
                //lmsAdapter(st3, 0.01, 1024, 0.5),
                //makeRLSAdapter(st3, 1, 1E-7),
                //st1,
                //lmsAdapter(st1, 0.01, 1024, 0.5),
                //makeRLSAdapter(st1, 1, 1E-7),
                );
    }
}

