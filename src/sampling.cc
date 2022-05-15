#include "Python/Python.h"
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "libs/pcg-cpp/include/pcg_random.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <random>

using namespace itensor;
using namespace std;
using std::vector;

int valueAtbit(int num, int bit)
{
    return (num >> bit) & 1;
}

vector<int> sample(MPS &psi)
{
    auto sites = siteInds(psi);
    auto N = length(psi);

    auto cps = vector<int>(N + 1);

    psi.position(1);
    for (int j = 1; j <= N; ++j)
    {
        Index sj = sites(j);

        auto PUp = ITensor(sj, prime(sj));

        PUp.set(1, 1, 1.0);

        Real prob_up = real(eltC(dag(prime(psi(j), "Site")) * PUp * psi(j)));

        int st = 0;
        if (Global::random() > prob_up)
            st = 1;
        cps.at(j) = st;

        auto upState = ITensor(sj);
        auto downState = ITensor(sj);

        upState.set(1, 1.0);
        downState.set(2, 1.0);

        // Project state
        ITensor jstate = (st == 1) ? downState : upState;
        if (j < N)
        {
            auto newA = psi(j + 1) * (dag(jstate) * psi(j));
            newA /= norm(newA);
            psi.set(j + 1, newA);
        }
        // Set site j tensor
        psi.set(j, jstate);
    }

    return cps;
}

int mwpm_python(int dx, int dy, int i_psi)
{
    int res;
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue, *sys, *path;

    // if (Stop == 1)
    // {
    //     Py_Initialize();
    // }

    sys = PyImport_ImportModule("sys");
    path = PyObject_GetAttrString(sys, "path");
    PyList_Append(path, PyUnicode_FromString("py_decoder/"));

    /* import */
    pName = PyUnicode_FromString("tcmain");
    pModule = PyImport_Import(pName);

    if (!pModule)
    {
        PyErr_Print();
        printf("ERROR in pModule\n");
        exit(1);
    }

    /* call python function */
    pFunc = PyObject_GetAttrString(pModule, "main_fun");

    /* build args */
    pArgs = PyTuple_New(3);
    PyTuple_SetItem(pArgs, 0, PyLong_FromLong(dx));
    PyTuple_SetItem(pArgs, 1, PyLong_FromLong(dy));
    PyTuple_SetItem(pArgs, 2, PyLong_FromLong(i_psi));

    /* call */
    pValue = PyObject_CallObject(pFunc, pArgs);

    res = PyLong_AsLong(pValue);

    // if (Stop == 3)
    // {
    //     Py_Finalize();
    // }

    return res;
}

// vector<int> numtovec(int &num)
// {
//     string binary = bitset<20>(num).to_string();
// }

double sample_measure(MPS &psi)
{
    auto N = length(psi);
    vector<int> vec = sample(psi);
    vector<int> vec_relabel(0);
    int num_array[] = {2, 1, 3, 0};
    for (int n : num_array)
    {
        for (int i = 0; i < N; i += 4)
        {
            vec_relabel.push_back(vec[n + i]);
        }
    }

    // int L = vec_relabel.size();
    int num_vec = 0;
    for (int i = 0; i < N; i++)
    {
        num_vec += vec_relabel[i] * 2 ^ (N - 1 - i);
    }
    int dx = N / 4;
    int dy = 2;
    int num_corr = mwpm_python(dx, dy, num_vec);

    int measure_lgz = 0;

    for (int i = 0; i < N / 4; i++)
    {
        int bit_i = N - 1 - i;
        measure_lgz += valueAtbit(num_corr, bit_i);
    }

    double res;

    if (measure_lgz % 2 == 0)
    {
        res = 1.0;
    }
    else
    {
        res = -1.0;
    }

    return res;
}

double measure_average(MPS &psi, int N_sample)
{
    double measure_a = 0;
    for (int i = 0; i <= N_sample; i++)
    {
        measure_a += sample_measure(psi) / N_sample;
    }
    return measure_a;
}
