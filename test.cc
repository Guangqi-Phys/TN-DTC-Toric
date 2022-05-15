#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
using namespace std;
using std::vector;

int main()
{
     // int num = 5;
     // int bit = 1;
     // int y = (num >> bit) & 1;
     // cout << y << endl;
     vector<int> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

     for (int i = 0; i < 12; i++)
          cout << vec[i] << endl;

     vector<int> vec_relabel(0);
     int num_array[] = {2, 1, 3, 0};
     for (int n : num_array)
     {
          // int n = num_array[i];
          // cout << "n = " << n << endl;
          for (int j = 0; j < 12; j += 4)
          {
               vec_relabel.push_back(vec[n + j]);
          }
     }

     for (int i = 0; i <= 12; i++)
          cout << i << endl;
     int L = vec_relabel.size();
     // cout << L << endl;

     return 0;
}
