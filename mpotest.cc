#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "libs/pcg-cpp/include/pcg_random.hpp"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

using namespace itensor;
using namespace std;
using std::vector;

int main(int argc, char *argv[])
{

    int N = 8;
    double theta0; // error rate
    // double theta0_prob;
    // double theta1; // stabilizer operator
    double theta2; // logical x operator
    int const n_t = 50;
    int const n_simu = 100;
    double measure0;
    double measure1;
    double measure2;
    auto sites = SpinHalf(N, {"ConserveQNs=", false});
    auto state = InitState(sites);
    double measur1list[n_t] = {0};
    double measur2list[n_t] = {0};
    string filename;

    // PrintData(state);
    // theta0 = 0.02 * M_PI;
    // theta1 = 0.10 * M_PI;
    theta2 = 0.5 * M_PI;

    pcg_extras::seed_seq_from<random_device> seed_source;

    pcg32 engine0(seed_source);
    pcg32 engine1(seed_source);

    uniform_real_distribution<double> dist0(-0.02, 0.02);

    // define the logical z operator
    auto logzmpo = AutoMPO(sites);
    logzmpo += pow(2, (N / 4)), "Sz", 3, "Sz", 7;
    auto Hlogz = toMPO(logzmpo);

    auto args = Args("Cutoff=", 1E-16, "MaxDim=", 200);

    // Make initial MPS psi to be in the all up state
    for (auto j : range1(N))
    {
        state.set(j, "Up");
    }

    auto psi = MPS(state);

    // applying stabilizer operators
    // for (int b = 0; b < N; b++)
    // {
    //     auto stbmpo = AutoMPO(sites);
    //     if (b % 4 == 0)
    //     {
    //         stbmpo += 4 * sin(theta1) * Cplx_i, "Sz", ((b % 12) + 1), "Sz", (((b + 2) % 12) + 1), "Sz", (((b + 3) % 12) + 1), "Sz", (((b + 4) % 12) + 1);
    //     }
    //     if (b % 4 == 1)
    //     {
    //         stbmpo += 4 * sin(theta1) * Cplx_i, "Sz", ((b % 12) + 1), "Sz", (((b + 1) % 12) + 1), "Sz", (((b + 2) % 12) + 1), "Sz", (((b + 4) % 12) + 1);
    //     }
    //     if (b % 4 == 2)
    //     {
    //         stbmpo += 4 * sin(theta1) * Cplx_i, "Sx", ((b % 12) + 1), "Sx", (((b + 2) % 12) + 1), "Sx", (((b + 3) % 12) + 1), "Sx", (((b + 4) % 12) + 1);
    //     }
    //     if (b % 4 == 3)
    //     {
    //         stbmpo += 4 * sin(theta1) * Cplx_i, "Sx", ((b % 12) + 1), "Sx", (((b + 1) % 12) + 1), "Sx", (((b + 2) % 12) + 1), "Sx", (((b + 4) % 12) + 1);
    //     }
    //     // note that need to do noPrime operation after applyMPO() function.
    //     auto Hstb = toMPO(stbmpo);
    //     psi = sum(cos(theta1) * psi, applyMPO(Hstb, psi, args).noPrime("Site"), args);
    //     psi.normalize();
    // }

    PrintData(psi);

    // apply logical x operator to state psi

    // auto logxmpo = AutoMPO(sites);
    // logxmpo += 4, "Sx", 7, "Sx", 8;
    // auto Hlogx = toMPO(logxmpo);
    // auto expHlogx = toExpH(logxmpo, theta2 * Cplx_i);
    // psi = applyMPO(expHlogx, psi, args).noPrime("Site");
    // psi.normalize();

    auto logxmpo = AutoMPO(sites);
    // auto logxmpo_t = AutoMPO(sites);
    logxmpo += 2 * sin(theta2) * Cplx_i, "Sx", 2;
    // // logxmpo_t += 0.3, "Id", 1, "Id", 2, "Id", 3, "Id", 4, "Id", 5, "Id", 6, "Id", 7, "Id", 8;
    auto Hlogx = toMPO(logxmpo);
    // auto Hlogx_t = toMPO(logxmpo_t);
    // psi = applyMPO(Hlogx_t, psi, args).noPrime("Site");
    // psi = 0.3 * psi;
    // psi = sum(applyMPO(Hlogx_t, psi, args).noPrime("Site"), applyMPO(Hlogx, psi, args).noPrime("Site"), args);
    psi = sum(cos(theta2) * psi, applyMPO(Hlogx, psi, args).noPrime("Site"), args);
    // psi.normalize();

    // for (int i = 7; i < 9; i++)
    // {
    //     auto logxmpo = AutoMPO(sites);
    //     logxmpo += sin(theta2) * Cplx_i, "Sx", i;
    //     auto Hlogx = toMPO(logxmpo);
    //     psi = sum(cos(theta2) * psi, applyMPO(Hlogx, psi, args).noPrime("Site"), args);
    //     psi.normalize();
    //     // auto expHlogx = toExpH(logxmpo, theta2 * Cplx_i);
    //     // psi = applyMPO(expHlogx, psi, args).noPrime("Site");
    // }

    PrintData(psi);

    return 0;
}
