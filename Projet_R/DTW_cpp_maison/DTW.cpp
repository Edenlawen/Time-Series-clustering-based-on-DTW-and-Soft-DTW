#include "include/DTW.h"
#include <bits/stdc++.h>

using namespace std;

void DTW(vector<double> P, vector<double> Q)
{
    // Initialisation de la matrice de taille x = P.size() et y = Q.size()
    vector<vector<double>> matrice(P.size(), vector<double>(Q.size()));

    // On calcul la valeur aux coordonnées [0][0] pour avoir une 1er valeur et pouvoir lancer le programme ensuite
    matrice[0][0] = abs(P[0] - Q[0]);

    /*
        Etape 1:
        On va mettre en place deux boucle for:
        - La première va permettre de remplir la 1er ligne à partir de l'index 1
        - La deuxième va permettre de remplir la 1er colonne à partir de l'index 1
    */
    cout << "Etape 1 qui commence" << endl;
    // Première Boucle
    for (int x = 1; x < Q.size(); x++)
    {
        matrice[0][x] = abs(P[0] - Q[x]) + matrice[0][x - 1];
    }
    // Deuxième Boucle
    for (int y = 1; y < P.size(); y++)
    {
        matrice[y][0] = abs(P[y] - Q[0]) + matrice[y - 1][0];
    }
    cout << "Etape 1 fini" << endl;

    /*
        Etape 2:
        On va calculer ligner par ligne à partir de i et j égal à 1.
        On va utiliser la formule ci-dessous pour calculer case par case.
    */
    cout << "Etape 2 qui commence" << endl;
    for (int i = 1; i < P.size(); i++)
    {
        for (int j = 1; j < Q.size(); j++)
        {
            matrice[i][j] = abs(P[i] - Q[j]) + min({matrice[i - 1][j], 2 * matrice[i - 1][j - 1], matrice[i][j - 1]});
        }
    }
    cout << "Etape 2 fini" << endl;

    /*
        Etape 3:
        Le but va être de trouver le chemin reliant la fin au début
    */
    cout << "Etape 3 qui commence" << endl;
    int coord_x = P.size() - 1;
    int coord_y = Q.size() - 1;

    vector<double> chemin_val;
    vector<string> chemin;
    string str = "";

    chemin_val.push_back(matrice[coord_x][coord_y]);
    str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
    chemin.push_back(str);

    bool cond = true;

    while (cond)
    {
        str = "";
        if (coord_x != 0 && coord_y != 0)
        {
            if ((matrice[coord_x - 1][coord_y] >= matrice[coord_x - 1][coord_y - 1]) && (matrice[coord_x][coord_y - 1] >= matrice[coord_x - 1][coord_y - 1]))
            {
                coord_x--;
                coord_y--;
                chemin_val.push_back(matrice[coord_x][coord_y]);
                str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
                chemin.push_back(str);
            }
            else
            {
                if (matrice[coord_x - 1][coord_y] >= matrice[coord_x][coord_y - 1])
                {
                    coord_y--;
                    chemin_val.push_back(matrice[coord_x][coord_y]);
                    str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
                    chemin.push_back(str);
                }
                else
                {
                    coord_x--;
                    chemin_val.push_back(matrice[coord_x][coord_y]);
                    str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
                    chemin.push_back(str);
                }
            }
        }
        else
        {
            if (coord_x == 0)
            {
                coord_y--;
                chemin_val.push_back(matrice[coord_x][coord_y]);
                str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
                chemin.push_back(str);
            }
            else
            {
                coord_x--;
                chemin_val.push_back(matrice[coord_x][coord_y]);
                str += "M[" + to_string(coord_x) + "][" + to_string(coord_y) + "]";
                chemin.push_back(str);
            }
        }
        if (coord_x == 0 && coord_y == 0)
            cond = false;
    }
    cout << "Etape 3 fini" << endl;

    /*
        Etape 4:
        Calcul de la distance
    */
    cout << "Etape 4 qui commence" << endl;
    double di_tot = 0;
    for (int i = 0; i < chemin_val.size(); i++)
        di_tot += chemin_val[i];

    double resultat = di_tot / chemin_val.size();
    cout << "Etape 4 fini\n" << endl;

    cout << "Le cout est de " << matrice[matrice.size()-1][matrice.size()-1] << "\n" << endl;

    // Chemin emprunté
    cout <<  "Position|Valeur"<< endl;
    for (int i = 0; i < chemin_val.size(); i++){
        cout << chemin[i] << "|" << chemin_val[i] <<endl;
    }

    cout << "\nLe distance est de " << resultat << "\n" << endl;

    // Code pour afficher la matrice
    for (int i = 0; i < matrice.size(); i++)
    {
        for (int j = 0; j < matrice[i].size(); j++)
        {
            cout << setw(4) << matrice[i][j] << " ";
        }
        cout << endl;
    }
}