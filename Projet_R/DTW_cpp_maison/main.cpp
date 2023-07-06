#include "include/DTW.h"
#include <bits/stdc++.h>

using namespace std;

int main()
{
    // Initialisation des jeu de données
    
    // vector<double> P = {1, 4, 5, 10, 9, 3};
    // vector<double> P = {1, 7, 3, 4, 1, 10, 5};
    // vector<double> P = {2, 1, 4, 5, 1,
                        // 7, 3, 4, 1, 10, 5
                        // };
    // vector<double> Q = {1, 7, 3, 4, 1, 10, 5};

    vector<double> P = {1,10,16,5};
    vector<double> Q = {1,4,13,9};

    // vector<double> test = P;

    // for (int i = 0; i < test.size(); i++)
    //     test[i] *= 2;

    cout << P.size()*Q.size()<<endl;

    DTW(P, Q);
    return 0;
}
