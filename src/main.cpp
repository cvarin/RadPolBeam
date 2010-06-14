#include <cmath>
#include <complex> 
#include <cstdio>
#include <iostream>
using namespace std;
typedef complex<double> dbcomplex;

// Bessel complexe : http://pagesperso-orange.fr/jean-pierre.moreau/c_bessel.html

/******************************************************************************/
/****************** Fonctions locales (prototypes) ****************************/
/******************************************************************************/
void print_stars();
void print_start_message();

/******************************************************************************/
/****************** Main ******************************************************/
/******************************************************************************/
int main(void)
{    
     /************* Message de départ *****************************************/
     print_start_message();
     
     /*************************************************************************/
//      int n = 10;
     dcomplex a(1.0,3.0);
     
     cout << spherical_j0(a) << endl;
     cout << spherical_j1(a) << endl;
     cout << spherical_j2(a) << endl;
     
     /************* Ligne d'étoiles *******************************************/
     print_stars();
     printf("\n");
}

/******************************************************************************/
/****************** Fonctions locales *****************************************/
/******************************************************************************/
void print_stars()
{printf("******************************************************************\n");}

/******************************************************************************/
void print_start_message()
{
     printf("\n");
     print_stars();
     printf("Programme pour tracer le champ d'un faisceau polarisé radialement.\n");
     printf("Ref: A. April, Opt. Lett. 331563 (2008).\n");
     printf("\n");
}

/****************** End of file ***********************************************/