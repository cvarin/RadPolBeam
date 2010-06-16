#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Ref pour fonctions de bessel sphériques : Arfken et Weber, "Mathematical methods
// for physicists", 4ième édition, section 11.7.

/******************************************************************************/
/****************** Fonctions locales (prototypes) ****************************/
/******************************************************************************/
complex double spherical_j0(complex double z);
complex double spherical_j1(complex double z);
complex double spherical_j2(complex double z);
void print_stars();
void print_start_message();

/******************************************************************************/
/****************** Main ******************************************************/
/******************************************************************************/
int main(int argc, char *argv[])
{    
     /************* Message de départ *****************************************/
     print_start_message();
     
     /*************************************************************************/
     const double pi = 3.14159265358979323846;
     const double k = 2.0*pi;
     const double zo = atof(argv[1])/k;
     const double wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+pow(k*zo,2))-1.0);
     
     /************* Variable radiale ******************************************/
     const int n = 1000;
     double rmax = 10.0;
     if(k*zo <= 100) rmax = 6.0;
     if(k*zo <= 70) rmax = 4.0;
     if(k*zo <= 30) rmax = 3.0;
     if(k*zo <= 10) rmax = 2.0;
     const double rmin = -rmax;
     double dr = (rmax - rmin)/(n-1);
     double r[n];
     r[0] = rmin;
     for(int i = 1; i<n; i++) r[i] = r[i-1] + dr;
     
     /************* Variables complexes ***************************************/
     complex double Z = 0.0 + zo*I;
     complex double kR[n], theta[n];
     for(int i = 0; i<n; i++)
     {
          kR[i] = k*csqrt(r[i]*r[i] + Z*Z);
          theta[i] = cacos((k*Z)/kR[i]);
     }
     
     /************* Champ électrique ******************************************/
     complex double Er[n],Ez[n];
     for(int i = 0; i<n; i++)
     {
          Er[i] = 0.5*spherical_j2(kR[i])*csin(2.0*theta[i]);
          Ez[i] = 2.0/3.0*(spherical_j0(kR[i]) + spherical_j2(kR[i])
                   *0.25*(1.0+3.0*ccos(2.0*theta[i])) );
     }
     
     /************* Intensité *************************************************/
     double wr[n],wz[n];
     for(int i = 0; i<n; i++)
     {
          wr[i] = cabs(Er[i])*cabs(Er[i]);
          wz[i] = cabs(Ez[i])*cabs(Ez[i]);
     }

     /************* Écrire dans un fichier ************************************/
     char nom[100];
     sprintf(nom,"intensite_kzo_%.3f.dat",(k*zo));
     FILE* file = fopen(nom,"w");
     
     for(int i = 0; i<n; i++)
          fprintf(file,"%f %f %f %f\n",r[i],wr[i],wz[i],wr[i]+wz[i]);
     fprintf(file,"\n");
     
     fclose(file);
     
     /************* Sortie à l'écran ******************************************/
     printf("k*zo = %f\n",zo*k);
     printf("wo = %f lambda\n",wo);
     printf("Sortie écrite dans %s\n",nom);
     printf("\n");
     
     /************* Ligne d'étoiles *******************************************/
     print_stars();
     printf("\n");
}

/******************************************************************************/
/****************** Fonctions locales *****************************************/
/******************************************************************************/
complex double spherical_j0(complex double z)
{
     return csin(z)/z;
}

/******************************************************************************/
complex double spherical_j1(complex double z)
{
     return csin(z)/(z*z) - ccos(z)/z;
}

/******************************************************************************/
complex double spherical_j2(complex double z)
{
     return (3.0/(z*z*z) - 1.0/z)*csin(z) - 3.0*ccos(z)/(z*z);
}

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