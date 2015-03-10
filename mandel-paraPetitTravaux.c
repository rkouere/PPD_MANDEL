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


RENDU
comment on a fait pour faire maitre esclave
Campagne de test (pas trop poussé mais il faut quand me la faire. Il faut ausis analyser ce que l'on a fait
Pour mesuerer, dans ping pong on a ce qu'il faut
Mesurer sur 2/4 machine et des images de taille differente.

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

#define X_SIZE  8192		/* dimension image */
#define Y_SIZE  6144
#define FILENAME "mandel.ppm"	/* image resultat */
#define TAG_LIMIT 1
#define TAG_RESULT 2
typedef struct {  
    int x_size, y_size;		/* dimensions */
    char *pixels;		/* matrice linearisee de pixels */
} picture_t;

/*------------------------------
  Les processus
  ------------------------------------------------------------*/
int self;			/* mon rang parmi les processus */
int procs;			/* nombre de processus */
#define MASTER 0

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

    /* ajout petits travaux */
    int nombre_de_bandes_total = 120; /* nb bandes que l'on va traiter */
    int nombreDeBandesTraitees = 0, i;
    int bande_en_cour_sur_proc[procs]; /* utilise pour pouvoir avoir un lien entre un processeur et la bande qu'il est en train de traiter (on ne pourrait sinon pas savoir quelle partie de l'image on doit mettre à jour lorsque le master reçoit des données d'un processeur). Commence à zero */
    double y_min_tmp, y_max_tmp;
    MPI_Request request[procs];
    int index;
    
    double startwtime = 0.0, endwtime;

    /* printf("start\n"); */


    com = MPI_COMM_WORLD;



    /* printf("MPI_COMM\n"); */
    MPI_Init (&argc, &argv);

    /* if(self == MASTER) */
    /*   startwtime = MPI_Wtime(); */

    /* printf("MPI_Init\n"); */
    MPI_Comm_size (com, &procs);
    MPI_Comm_rank (com, &self);
    /* printf("procs %d\n", procs); */
    /* printf("self %d\n", self); */

    parse_argv(argc, argv,
    	       &n_iter,
    	       &x_min, &x_max, &y_min, &y_max,
    	       &x_size, &y_size, &pathname);


    /* on calcul le pas de chaque bandes */
    pas_y = (y_max - y_min)/nombre_de_bandes_total;
    /* on initialise une structure de données qui permettra de stoquer chaque bande... c'est un peut comme un fichier temporaire */
    init_picture (&pict, x_size, y_size/nombre_de_bandes_total);
    if(self == MASTER) {
      /* initialisation de la structure de donné utilisé pour stocquer les informaitons de l'image finale */
      init_picture (&pictFinal, x_size, y_size);
      /* phase d'initialisaiton */
      /* pour chaque processeur utilisable : on envoit le y_min necessaire au calcul et on incremente le nombre de bande traite */
      for(i = 1; i < procs; i ++) {
	/* on calcul le y_min que l'on va envoyer en fonction du nombre de bande déjà traité */
	y_min_tmp = y_min + pas_y*(nombre_de_bandes_total-nombreDeBandesTraitees-1);
	/* on fait le lien numero de bande > processeur pour la reception */
	bande_en_cour_sur_proc[i] = nombreDeBandesTraitees;
	nombreDeBandesTraitees++;
	/* on envoit la données à chaque processeur */
	MPI_Send(&y_min_tmp, 1, MPI_DOUBLE, i, TAG_LIMIT, com);
      }
      /* notre compteur est egale au nombre de procs */
      
      /* gestion de la charge */
      /* tant qu'il y a des bandes à traiter */
      while(nombreDeBandesTraitees < nombre_de_bandes_total) {
	/* on recupere les donnees de la bande traite par le processeur */
	MPI_Recv(pict.pixels, x_size*(y_size/nombre_de_bandes_total), MPI_CHAR, MPI_ANY_SOURCE, TAG_RESULT, com, &status);
	/* on copie ces données à l'endroit ou il faut dans la structure contenant les données de l'image final */
	memcpy(pictFinal.pixels+(bande_en_cour_sur_proc[status.MPI_SOURCE])*x_size*(y_size/nombre_de_bandes_total), pict.pixels, (x_size*(y_size/nombre_de_bandes_total)));
	/* on calcul le y_min que l'on va envoyer en fonction du nombre de bande déjà traité */
	y_min_tmp = y_min + pas_y*(nombre_de_bandes_total-nombreDeBandesTraitees-1);
	/* on fait le lien numero de bande > processeur pour la reception */
	bande_en_cour_sur_proc[status.MPI_SOURCE] = nombreDeBandesTraitees;
	nombreDeBandesTraitees++;
	/* on envoit les nouvelles données au processeur libre */
	MPI_Send(&y_min_tmp, 1, MPI_DOUBLE, status.MPI_SOURCE, TAG_LIMIT, com);
      }

      /* il reste encore à recevoir les dernièrs calculs */
      for(i = 1; i < procs; i++) {
      	MPI_Recv(pict.pixels, x_size*(y_size/nombre_de_bandes_total), MPI_CHAR, MPI_ANY_SOURCE, TAG_RESULT, com, &status);
      	memcpy(pictFinal.pixels+(bande_en_cour_sur_proc[status.MPI_SOURCE])*x_size*(y_size/nombre_de_bandes_total), pict.pixels, (x_size*(y_size/nombre_de_bandes_total)));
      }
      y_min_tmp = -1;

      /* on utilise cette technique un peut "crade" pour indiquer aux processeur qu'ils n'ont plus de travail à réaliser... il aurait ete plus propre d'utiliser un autre tag pour cela et faire un broadcast */
      for(i = 1; i < procs; i++) {
	MPI_Send(&y_min_tmp, 1, MPI_DOUBLE, i, TAG_LIMIT, com);
      }
      
    }

    /* calculs des bandes */
    if(self != MASTER) {
      /* tant qu'il y a des bandes à traiter (tant que le y_min != -1) */
      while(1) {
	MPI_Recv(&y_min_tmp, 1, MPI_DOUBLE, MASTER, TAG_LIMIT, com, &status);
	if(y_min_tmp == -1)
	  break;
	/* le processeur calcul sa bande et la renvoit au master */
	compute (&pict, n_iter, x_min, x_max, y_min_tmp, y_min_tmp + pas_y); 
	MPI_Send(pict.pixels, x_size*(y_size/nombre_de_bandes_total), MPI_CHAR, MASTER, TAG_RESULT, com);
      }
      
    }

    
    sprintf(localpathName, "yo%d.ppm", self);

    if(self == MASTER) {
      /* printf("saving picture %s\n", localpathName); */
      save_picture (&pictFinal, localpathName);
    }
    
    /* if(self == MASTER) { */
    /*   endwtime = MPI_Wtime(); */
    /*   printf("time = %f\n", endwtime-startwtime);	        */
    /* } */
    /* on ferme ! */
    MPI_Finalize ();
    exit(EXIT_SUCCESS);
}
