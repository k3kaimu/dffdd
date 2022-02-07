module dffdd.detector.primitives;


interface IDetector(X, Y)
{
    alias InputElementType = X;
    alias OutputElementType = Y;

    size_t inputLength() const;
    size_t outputLength() const;
    Y[] detect(in X[], return ref Y[]);
}
