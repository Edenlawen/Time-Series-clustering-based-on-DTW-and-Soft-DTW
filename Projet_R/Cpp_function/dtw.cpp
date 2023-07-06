#include <Rcpp.h>
#include <bits/stdc++.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix costMatrix(vector<double> P, vector<double> Q, bool verbose = false){
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
  if(verbose==true){cout << "Etape 1 qui commence" << endl;}
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
  if(verbose==true){cout << "Etape 1 fini" << endl;}
  
  /*
   Etape 2:
   On va calculer ligner par ligne à partir de i et j égal à 1.
   On va utiliser la formule ci-dessous pour calculer case par case.
   */
  if(verbose==true){cout << "Etape 2 qui commence" << endl;}
  for (int i = 1; i < P.size(); i++)
  {
    for (int j = 1; j < Q.size(); j++)
    {
      int d = abs(P[i] - Q[j]);
      matrice[i][j] = min({d + matrice[i - 1][j], 2 * d + matrice[i - 1][j - 1], d + matrice[i][j - 1]});
      /* Méthode de la prof
      matrice[i][j] = d + min({matrice[i - 1][j], 2 * matrice[i - 1][j - 1], matrice[i][j - 1]});
      */
    }
  }
  if(verbose==true){cout << "Etape 2 fini" << endl;}
  
  NumericMatrix myMatrix(P.size(), Q.size());
  for(int i = 0; i<P.size(); i++){
    for(int j = 0; j < Q.size(); j++){
      myMatrix(i,j) = matrice[i][j];
    }
  }
  
  return myMatrix;
}

// [[Rcpp::export]]
double costFin(vector<double> P, vector<double> Q, bool verbose = false){
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
  if(verbose==true){cout << "Etape 1 qui commence" << endl;}
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
  if(verbose==true){cout << "Etape 1 fini" << endl;}
  
  /*
   Etape 2:
   On va calculer ligner par ligne à partir de i et j égal à 1.
   On va utiliser la formule ci-dessous pour calculer case par case.
   */
  if(verbose==true){cout << "Etape 2 qui commence" << endl;}
  for (int i = 1; i < P.size(); i++)
  {
    for (int j = 1; j < Q.size(); j++)
    {
      int d = abs(P[i] - Q[j]);
      matrice[i][j] = min({d + matrice[i - 1][j], 2 * d + matrice[i - 1][j - 1], d + matrice[i][j - 1]});
      /* Méthode de la prof
      matrice[i][j] = d + min({matrice[i - 1][j], 2 * matrice[i - 1][j - 1], matrice[i][j - 1]});
      */
    }
  }
  if(verbose==true){cout << "Etape 2 fini" << endl;}
  

  
  return matrice[P.size() - 1][Q.size() - 1];
}