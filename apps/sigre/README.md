# SigRe: Signal Regenerator

## Usage

```sh
$ sigre <nTrain> <nRegen> <Regenerator Type> [Regenerator Parameters...] @ <Basis Function Builders...>
```

## Examples

```sh
$ sigre 1000 80 PHammLS 32 @ OFDMUpSampler 64 16 8 % BFList X[1,0] L[0,1] X[3,0] M[0,1,*1] % OFDMDownSampler 64 16 8 % Orthonormalizer
```
