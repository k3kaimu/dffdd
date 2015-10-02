import std.stdio;

import dffdd.filter.ls;

void main()
{
    cfloat[6][5] mat;
    cfloat[6] vec;
    float[30] s;

    mat[0][0] = 0i -0.09; mat[1][0] = 0i +0.14; mat[2][0] = 0i -0.46; mat[3][0] = 0i +0.68; mat[4][0] = 0i +1.29;
    mat[0][1] = 0i -1.56; mat[1][1] = 0i +0.20; mat[2][1] = 0i +0.29; mat[3][1] = 0i +1.09; mat[4][1] = 0i +0.51;
    mat[0][2] = 0i -1.48; mat[1][2] = 0i -0.43; mat[2][2] = 0i +0.89; mat[3][2] = 0i -0.71; mat[4][2] = 0i -0.96;
    mat[0][3] = 0i -1.09; mat[1][3] = 0i +0.84; mat[2][3] = 0i +0.77; mat[3][3] = 0i +2.11; mat[4][3] = 0i -1.27;
    mat[0][4] = 0i +0.08; mat[1][4] = 0i +0.55; mat[2][4] = 0i -1.13; mat[3][4] = 0i +0.14; mat[4][4] = 0i +1.74;
    mat[0][5] = 0i -1.59; mat[1][5] = 0i -0.72; mat[2][5] = 0i +1.06; mat[3][5] = 0i +1.24; mat[4][5] = 0i +0.34;

    vec[0] = 0i +7.4;
    vec[1] = 0i +4.2;
    vec[2] = 0i -8.3;
    vec[3] = 0i +1.8;
    vec[4] = 0i +8.6;
    vec[5] = 0i +2.1;

    int rank, info;
    LAPACKE_cgelss(102, 6, 5, 1, mat[0].ptr, 6,
               vec.ptr, 6, s.ptr, 0.01f,
               &rank/*, work.ptr, 1024, rwork.ptr, &info*/).writeln;

    writeln(vec);
    writeln(s);
    writeln(rank);
    writeln(info);
}
