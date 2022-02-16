module dffdd.mod.filterbank;

import std.complex;
import dffdd.filterbank.dftfilterbank;

final class OversampledDFTFilterBankModulator(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(size_t nchannel, in C[] prototype)
    {
        this._analysis = new OversampledDFTAnalysisFilterBank!C(nchannel, prototype);
        this._synthesis = new OversampledDFTSynthesisFilterBank!C(nchannel, prototype);
    }


    size_t symInputLength() const @property
    {
        return this._analysis.numOfChannel();
    }


    size_t symOutputLength() const @property
    {
        return this._analysis.numOfChannel() / 2;
    }


    OutputElementType[] modulate(in InputElementType[] inputs, OutputElementType[] outputs)
    in {
        assert(inputs.length % this.symInputLength == 0);
        assert(outputs.length == (inputs.length / this.symInputLength * this.symOutputLength));
    }
    do {
        return this.modulate(inputs, outputs);
    }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in {
        assert(inputs.length % this.symInputLength == 0);
    }
    do {
        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / 2;
        
        foreach(i; 0 .. inputs.length / this.symInputLength) {
            _analysis(
                inputs[i * this.symInputLength .. (i+1) * this.symInputLength],
                outputs[i * this.symOutputLength .. (i+1) * this.symOutputLength]);
        }

        return outputs;
    }


    InputElementType[] demodulate(in OutputElementType[] inputs, InputElementType[] outputs)
    in {
        assert(inputs.length % this.symOutputLength == 0);
        assert(outputs.length == (inputs.length / this.symOutputLength * this.symInputLength));
    }
    do {
        return this.demodulate(inputs, outputs);
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in {
        assert(inputs.length % this.symOutputLength == 0);
    }
    do {
        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. inputs.length / this.symOutputLength) {
            _synthesis(
                inputs[i * this.symOutputLength .. (i+1) * this.symOutputLength],
                outputs[i * this.symInputLength .. (i+1) * this.symInputLength]
            );
        }

        return outputs;
    }


  private:
    OversampledDFTAnalysisFilterBank!C _analysis;
    OversampledDFTSynthesisFilterBank!C _synthesis;
}

unittest
{
    import std.algorithm : map;
    import std.range : array;
    import dffdd.filterbank.protofilter : designKaiserBesselDrived;
    import dffdd.window : kaiserBetaFromStopbandAttdB;


    alias R = float;
    alias C = Complex!R;
    immutable size_t nCH = 32,
                     nTap = 32;

    auto proto = designKaiserBesselDrived!R(nCH/2, nTap, kaiserBetaFromStopbandAttdB(200), 1).map!(a => complex!R(a, 0)).array();

    auto mod = new OversampledDFTFilterBankModulator!C(nCH, proto);
}

