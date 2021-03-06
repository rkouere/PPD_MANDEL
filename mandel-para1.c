/* ------------------------------
   $Id: mandel-seq.c,v 1.2 2008/03/04 09:52:55 marquet Exp $
   ------------------------------------------------------------
   
   Affichage de l'ensemble de Mandelbrot.
   Version sequentielle.
   
*/

/*
Etapes
1. Version parallele
Pb 2d : on va le decouper en 1 dimension.
-> decoupage statique (a priori) horizontal ou vertical (il y a une bonne solution et une moins bonne)

2. 
a) Découpage avec de petits travaux
-> on va faire des fines bandes vertical et on les envoit. 
+ meilleur equilibre de la charge
- il faut gerer le partage :
-- system type maitre esclave
-- autres solutions

b) Dans le maitre esclave + recouvrir les communications
M -> travail à faire (bout d'image) -> E
E -> envois le résultat -> M
-->> envoyer message j'ai fini et ensuite envoyer les résultats (ce qui permet au maître d'envoyer une nouvelle image en même temps qu'à l'escalve d'envoyer ses résultats)

2bis.
travaux de taille variable
On envoit au debut la moitié du travail
Puis on envoit la moitié du travail restant
Puis on envoit la moitier de la moitier etc...

2ter
Solution distribue




 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "mpi.h"
#include <string.h>

/* Valeur par defaut des parametres */
#define N_ITER  255		/* nombre d'iterations */

#define X_MIN   -1.78		/* ensemble de Mandelbrot */
#define X_MAX   0.78
#define Y_MIN   -0.961
#define Y_MAX   0.961

#define X_SIZE  2048		/* dimension image */
#define Y_SIZE  1536
#define FILENAME "mandel.ppm"	/* image resultat */

typedef struct {  
    int x_size, y_size;		/* dimensions */
    char *pixels;		/* matrice linearisee de pixels */
} picture_t;

/*------------------------------
  Les processus
  ------------------------------------------------------------*/
int self;			/* mon rang parmi les processus */
int procs;			/* nombre de processus */
#define PROC_NULL 0

static void 
usage() 
{
    fprintf(stderr, "usage : ./mandel [options]\n\n");
    fprintf(stderr, "Options \t Signification \t\t Val. defaut\n\n");
    fprintf(stderr, "-n \t\t Nbre iter. \t\t %d\n", N_ITER);
    fprintf(stderr, "-b \t\t Bornes \t\t %f %f %f %f\n",
	    X_MIN, X_MAX, Y_MIN, Y_MAX);
    fprintf(stderr, "-d \t\t Dimensions \t\t %d %d\n", X_SIZE, Y_SIZE);
    fprintf(stderr, "-f \t\t Fichier \t\t %s\n", FILENAME);

    exit(EXIT_FAILURE);
}

static void 
parse_argv (int argc, char *argv[], 
	    int *n_iter, 
	    double *x_min, double *x_max, double *y_min, double *y_max, 
	    int *x_size, int *y_size, 
	    char **path) 
{
    const char *opt = "b:d:n:f:";
    int c;
    
    /* Valeurs par defaut */
    *n_iter = N_ITER;
    *x_min  = X_MIN;
    *x_max  = X_MAX;
    *y_min  = Y_MIN;
    *y_max  = Y_MAX;
    *x_size = X_SIZE;
    *y_size = Y_SIZE;
    *path   = FILENAME;

    /* Analyse arguments */
    while ((c = getopt(argc, argv, opt)) != EOF) {    
	switch (c) {      
	    case 'b': 		/* domaine */
		sscanf(optarg, "%lf", x_min);
		sscanf(argv[optind++], "%lf", x_max);
		sscanf(argv[optind++], "%lf", y_min);
		sscanf(argv[optind++], "%lf", y_max);
		break;
	    case 'd':		/* largeur hauteur */
		sscanf(optarg, "%d", x_size);
		sscanf(argv[optind++], "%d", y_size);
		break;
	    case 'n':		/* nombre d'iterations */
		*n_iter = atoi(optarg);
		break;
	    case 'f':		/* fichier de sortie */
		*path = optarg;
		break;
	    default:
		usage();
	}
    }  
}

static void 
init_picture (picture_t *pict, int x_size, int y_size)
{  
    pict->y_size = y_size;
    pict->x_size = x_size;
    pict->pixels = malloc(y_size * x_size); /* allocation espace memoire */
} 

/* Enregistrement de l'image au format ASCII .ppm */
static void 
save_picture (const picture_t *pict, const char *pathname)
{  
    unsigned i;
    FILE *f = fopen(pathname, "w");  

    fprintf(f, "P6\n%d %d\n255\n", pict->x_size, pict->y_size); 
    for (i = 0 ; i < pict->x_size * pict->y_size; i++) {
	char c = pict->pixels[i];
	fprintf(f, "%c%c%c", c, c, c); /* monochrome blanc */
    }
    
    fclose (f);
}

static void
compute (picture_t *pict,
	 int nb_iter,
	 double x_min, double x_max, double y_min, double y_max)
{  
    int pos = 0;
    int iy, ix, i;
    double pasx = (x_max - x_min) / pict->x_size, /* discretisation */
	   pasy = (y_max - y_min) / pict->y_size; 

    /* Calcul en chaque point de l'image */
    for (iy = 0 ; iy < pict->y_size ; iy++) {
	for (ix = 0 ; ix < pict->x_size; ix++) {
	    double a = x_min + ix * pasx,
		b = y_max - iy * pasy,
		x = 0, y = 0;      
	    for (i = 0 ; i < nb_iter ; i++) {
		double tmp = x;
		x = x * x - y * y + a;
		y = 2 * tmp * y + b;
		if (x * x + y * y > 4) /* divergence ! */
		    break; 
	    }
      
	    pict->pixels[pos++] = (double) i / nb_iter * 255;    
	}
    }
}

int 
main (int argc, char *argv[]) 
{
  
    int n_iter,			/* degre de nettete  */  
      x_size, y_size;		/* & dimensions de l'image */  
    double x_min, x_max, y_min, y_max, pas_y; /* bornes de la representation */
    
    int segSize = 0; /* taille d'un segment */
    char *pathname, localpathName[100];		/* fichier destination */
    picture_t pict, pictFinal;

    MPI_Comm com;			/* un/le communicateur */
    MPI_Status status;		/* un status des receptions de message */
    printf("start\n");

    com = MPI_COMM_WORLD;
    printf("MPI_COMM\n");
    MPI_Init (&argc, &argv);
    printf("MPI_Init\n");
    MPI_Comm_size (com, &procs);
    MPI_Comm_rank (com, &self);
    printf("procs %d\n", procs);
    printf("self %d\n", self);

    parse_argv(argc, argv,
    	       &n_iter,
    	       &x_min, &x_max, &y_min, &y_max,
    	       &x_size, &y_size, &pathname);

    

    pas_y = (y_max - y_min)/procs;
    init_picture (&pict, x_size, y_size/procs);
    if(self == PROC_NULL)
      init_picture (&pictFinal, x_size, y_size);

    /* on fait procs-self-1 car on veut commence d'en bas */
    compute (&pict, n_iter, x_min, x_max, y_min + pas_y*(procs-self-1), y_min + pas_y*(procs-self));

    
    
    MPI_Gather(pict.pixels, x_size*(y_size/procs), MPI_CHAR, pictFinal.pixels, x_size*(y_size/procs), MPI_CHAR, 0, com); 
    
    sprintf(localpathName, "yo%d.ppm", self);

    if(self == PROC_NULL) {
      printf("saving picture %s\n", localpathName);
      save_picture (&pictFinal, localpathName);
    }
    
    /* on ferme ! */
    MPI_Finalize ();
    exit(EXIT_SUCCESS);
}
