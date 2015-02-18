/*- * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  $Id: df.c,v 1.2 2008/02/26 15:52:21 marquet Exp $
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  File:  df.c

  Utilisation : df filein fileout
  
  Difference finies 1-D en MPI

  La méthode des différences finies consiste à remplacer les dérivées
  du premier ou second ordre apparaissant dans un problème continu par
  des différences finies qui s'appuient sur un nombre fini de points
  discrets ou noeuds du maillage.

  Ce programme est une version preliminaire a modifier

  1. Utiliser des communications collectives
     a. pour la diffusion de la taille du tableau.
     b. pour la diffusion des valeurs du tableau.
     c. pour le rassemblement des résulats.
     d. pour le rassemblement de l'erreur.
  2. Recouvrir les communications entre voisins par du calcul.
     a. envoyer de manière asynchone les limites du domaine local.
     b. réaliser le calcul sur le centre du domaine local.
     c. recevoir les limites des domaines des voisins.
     d. terminer les calculs locaux. 
  3. Passer au 2-D
     a. modifier la fonction de calcul.  
     b. quel decoupage du domaine / distribution des données ?
     c. mise en oeuvre à l'aide de topologies de processus.
     d. mise en oeuvre à l'aide de types dérivés MPI.
  
  Ce programme est destiné à l'étude de MPI. On ne peut en faire une
  application performante sans des modifications majeures.
  
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static const char cvsid[] = "$Id: df.c,v 1.2 2008/02/26 15:52:21 marquet Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/param.h>
#include <sys/time.h>
#include <errno.h>
#include "mpi.h"

#define EPSILON 0.05		/* converge quand erreur < EPSILON */


/*------------------------------
  Petits outils
  ------------------------------------------------------------*/
static void
pwarn (const char *mess)
{
  printf ("Warning: %s\n", mess);
  fflush (stdout) ; 
}

static void
pfatal (const char *mess)
{
  perror (mess);
  MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
}

/*------------------------------
  Les processus
  ------------------------------------------------------------*/
int self;			/* mon rang parmi les processus */
int procs;			/* nombre de processus */
int rnbr, lnbr;			/* mes voisins de droite et de gauche */

#define PROC_NULL 0

/*------------------------------
  Des etiquettes pour les messages
  ------------------------------------------------------------*/
#define TAG_BCAST_SIZE  12	/* de PROC_NULL vers proc */
#define TAG_BCAST_VAL   13
#define TAG_RNEIGHBOR	14	/* vers le voisin de droite */
#define TAG_LNEIGHBOR	15	/* ou de gauche  */
#define TAG_COLLECT_VAL	16	/* de proc vers PROC_NULL */
#define TAG_COLLECT_ERR 17	/* collecte de l'erreur */

/*------------------------------
  Definition de tableaux
  ------------------------------------------------------------*/

typedef struct struct_array_t array_t;
struct struct_array_t {
  int size;			/* une taille */
  float *val;			/* et des valeurs */
};

static array_t
array_alloc (int size)
{
  array_t a;
  a.size = size;
  a.val = (float *) malloc (size * sizeof(float));
  if (! a.val) 
    pfatal ("Allocation tableau");
  return a;
}

static void
array_free (array_t *ap)
{
  ap->size = 0;
  free (ap->val);
}

static void
array_swap (array_t *ap1, array_t *ap2)
{
  float *pf;

  if (ap1->size != ap2->size) 
    pfatal ("Echange tableau");
  
  pf = ap1->val; ap1->val = ap2->val; ap2->val = pf;
  
}

/*------------------------------
  I/O de tableaux 
  ------------------------------------------------------------*/

/* un tableau est stocké dans un fichier sous la forme d'un entier (sa
   taille) et d'une suite de valeurs float. Une valeur par ligne. */
   
static array_t
array_read (const char *filename)
{
  FILE *fd;
  array_t a;
  int size;
  int i, s;
  
  fd = fopen (filename, "r");
  if (!fd)
    pfatal ("Lecture tableau, ouverture du fichier");
  
  s = fscanf (fd, "%d\n", &size);
  if (s !=1) 
    pfatal ("Lecture tableau, taille du tableau");

  a = array_alloc (size);
  
  for (i=0 ; i<a.size ; i++) {
    s = fscanf (fd, "%e\n", a.val+i);
    if (s != 1) 
      pfatal ("Lecture tableau, manque de valeur");
  }
  
  s = fclose (fd);
  if (s)
    pfatal ("Lecture tableau, fermeture fichier");

  return a;
}

static void
array_write (array_t a, const char *filename)
{
  int i, s;
  FILE *fd;
  
  fd = fopen (filename, "w");
  if (!fd) 
    pfatal ("Ecriture tableau, ouverture du fichier");

  fprintf (fd, "%d\n", a.size);

  for (i=0 ; i<a.size; i++) 
    fprintf (fd, "%e\n", a.val[i]);

  s = fclose (fd);
  if (s)
    pfatal ("Ecriture tableau, fermeture fichier");
}

static void
array_print (array_t a, const char *name)
{
#define VAL_PER_LINE 5
  int i;
  
  printf ("Tableau %s (%d valeurs) :\n", name, a.size);
  for (i=0 ; i<a.size; i++) {
    printf ("%3g ", a.val[i]);
    if (i%VAL_PER_LINE == VAL_PER_LINE-1) printf ("\n");
  }
  printf ("\n--\n");
  fflush (stdout); 
}

/*------------------------------
  Calcul sur les tableaux
  ------------------------------------------------------------*/

/* retourne la différence maximale en valeur absolue entre les valeurs
   des deux tableaux pour les [length] éléments à compter de l'élément
   d'indice [minloc].  */ 
static float
max_error (array_t a1, array_t a2, int minloc, int length)
{
  float max = 0.0;
  int i;

  for (i=minloc ; i<length+minloc ; i++) {
    float diff = fabs (a1.val[i]-a2.val[i]);
    if (diff > max) max = diff;
  }
  
  return max;
}

/* calcule dans [new] la nouvelle valeur de [old] pour [length]
   éléments à compter de l'élément d'indice [minloc]. accède donc en
   lecture aux éléments [minloc-1] à [length+minloc+1]. */
static float
compute (array_t old, array_t new, int minloc, int length)
{
  int i;

  for (i=minloc ; i<length+minloc ; i++) {
    new.val[i] = (2*old.val[i] + old.val[i-1] + old.val[i+1])/4;
  }
}

/*------------------------------
  Programme
  ------------------------------------------------------------*/
int
main (int argc, char *argv[])
{
  MPI_Comm com;			/* un/le communicateur */
  MPI_Status status;		/* un status des receptions de message */
  array_t a_io;			/* pour les io */
  array_t a_local;		/* la partie locale du tableau + 2 elements */
  array_t a_new;		/* et sa nouvelle valeur */
  int gsize, lsize;		/* nombre d'elements geres globalement
				   et localement */
  float gerror, lerror;		/* erreurs globale et locale */
  int iter=0;			/* nombre d'iterations avant convergenece */
  int p;
  char *filein, *fileout;	/* les fichiers d'entree et de sortie */
  
  /* initialisations MPI */
  com = MPI_COMM_WORLD;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (com, &procs);
  MPI_Comm_rank (com, &self);

  /* les arguments de la commande */
  if (argc != 3) {
    if (self == PROC_NULL) 
      printf ("Usage : %s filein fileout\n", argv[0]); 
    MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if (self == PROC_NULL) {
    filein=argv[1] ; fileout=argv[2]; 
  }
	  
  /* lecture d'un tableau par PROC_NULL */
  if (self == PROC_NULL) {
    a_io = array_read (filein);
    gsize = a_io.size;
  }
  
  /* taille du tableau local */
  if (self == PROC_NULL) {
    /* taille multiple de procs */
    if (gsize % procs) {
      pwarn ("Taille du tableau incompatible, tronquée");
      gsize = procs * (gsize / procs);
    }
    lsize = gsize/procs;
  }
  
  /* diffusion de la taille de ce tableau aux autres */
  if (self == PROC_NULL) {
    for (p=1 ; p<procs ; p++)
      MPI_Send (&lsize, 1, MPI_INT, p, TAG_BCAST_SIZE, com);
  } else {
    MPI_Recv (&lsize, 1, MPI_INT, PROC_NULL, TAG_BCAST_SIZE, com, &status);
  }

  /* allocation pour le tableau local. l'élément d'indice [0] est la
     valeur qui sera reçue du voisin de gauche. l'élément d'indice
     [lsize+1] est la valeur qui sera reçue du voisin de droite. */
  a_local = array_alloc (lsize+2);
  a_new = array_alloc (lsize+2);
  
  /* diffusion des valeurs du tableau */
  if (self == PROC_NULL) {
    for (p=1 ; p<procs ; p++)
      MPI_Send (a_io.val + (p*lsize),
		lsize, MPI_FLOAT,
		p, TAG_BCAST_VAL, com);
    /* diffusion de PROC_NULL a PROC_NULL : recopie locale du tableau */
    memcpy (a_local.val+1, a_io.val, lsize * sizeof (float)); 
  } else {
    MPI_Recv (a_local.val+1, lsize, MPI_FLOAT,
	      PROC_NULL, TAG_BCAST_VAL, com, &status);
  }

  /* mes voisins */
  lnbr = (self+procs-1) % procs ;
  rnbr = (self+1) % procs;
  
  /* boucle de calcul jusqu'à convergence */
  do {
    iter++;
    
    /* échange de valeurs avec les voisins */
    MPI_Send (a_local.val+1, 1, MPI_FLOAT,
	      lnbr, TAG_LNEIGHBOR, com);
    MPI_Recv (a_local.val+lsize+1, 1, MPI_FLOAT,
	      rnbr, TAG_LNEIGHBOR, com, &status);
    MPI_Send (a_local.val+lsize, 1, MPI_FLOAT,
	      rnbr, TAG_RNEIGHBOR, com);
    MPI_Recv (a_local.val, 1, MPI_FLOAT,
	      lnbr, TAG_RNEIGHBOR, com, &status);

    /* calcul de la nouvelle génération */
    compute (a_local, a_new, 1, lsize);
    
    /* calcul de l'erreur locale */
    lerror = max_error (a_local, a_new, 1, lsize);

    /* remplacement de génération */
    array_swap (&a_local, &a_new);
    
    /* calcul de l'erreur globale */
    MPI_Allreduce (&lerror, &gerror, 1, MPI_FLOAT, MPI_MAX, com);
    
  } while (gerror > EPSILON) ;
  
  /* centralisation des tableaux par le processus PROC_NULL */
  if (self == PROC_NULL) {
    for (p=1 ; p<procs ; p++)
      MPI_Recv (a_io.val + (p*lsize),
		lsize, MPI_FLOAT,
		p, TAG_COLLECT_VAL, com, &status);
    /* diffusion de PROC_NULL a PROC_NULL : recopie locale du resultat */
    memcpy(a_io.val, a_local.val+1, lsize * sizeof (float)); 
  } else {
    MPI_Send (a_local.val+1, lsize, MPI_FLOAT,
	      PROC_NULL, TAG_COLLECT_VAL, com);
  }

  /* sauvegarde du resultat par le processus PROC_NULL */
  if (self == PROC_NULL) {
    array_write (a_io, fileout) ;
    printf ("Convergence en %d iterations\n", iter); 
    printf ("Resultat sauvegarde dans %s\n", fileout) ; 
  }
  
  /* on ferme ! */
  MPI_Finalize ();
  exit (EXIT_SUCCESS) ; 
}
