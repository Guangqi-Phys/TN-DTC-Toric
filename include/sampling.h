#ifndef SAMPLING_HPP_
#define SAMPLING_HPP_

using namespace itensor;
using namespace std;

int valueAtbit(int num, int bit);
vector<int> sample(MPS &psi);
double sample_measure(MPS &psi);
double measure_average(MPS &psi, int N_sample);
int mwpm_python(int dx, int dy, int i_psi);

#endif // SAMPLING_HPP_
