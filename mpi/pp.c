/*- * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  $Id: pp.c,v 1.1 2008/02/26 14:27:53 marquet Exp $
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  File:  pp.c
  
  mesure des performances des communications sur le schema de
  communication ping pong.  

  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static const char cvsid[] = "$Id: pp.c,v 1.1 2008/02/26 14:27:53 marquet Exp $" ;

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <values.h>
#include <sys/param.h>
#include <sys/time.h>
#include "mpi.h"

#define ITER 50
#define SIZE 125000

int
main (int argc, char *argv[]) 
{
  int i,
    self,			/* moi */
    procs,     			/* nb processus */
    bytes; 			/* taille messages */
  int
    dest=1,			/* processus destinataire */
    root=0;			/* processus emmeteur */
  int
    tag=52,			/* etiquette des messages */
    tagack=53; 			/* etiquette des messages retournes */
  double *array;		/* buffer du message */
  double			/* timers */
    start, end,
    elapsed,
    min=MAXDOUBLE,
    max=0.0, 
    dummy, overhead,
    total=0.0;
  char hostname [MAXHOSTNAMELEN]; /* nom du noeud */
  MPI_Status status;
  MPI_Request request;

  /* construction du buffer */
  bytes = sizeof (double) * SIZE;
  array = (double *) malloc (bytes);
  for (i=0; i<SIZE; i++) array[i] = 0.0;

  /* initialisation MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &self);

  /* sur quel noeud sommes nous */
  gethostname (hostname, MAXHOSTNAMELEN);

  /* coucou */
  printf( "Processus %d sur le noeud %s\n", self, hostname);

  /* info */ 
  if (self == root) printf("%d iterations, %d bytes\n", ITER, bytes);


  for (i =0 ; i<ITER; i++) {
    if (self == root) {
      /* pour test integrite */
      array [SIZE-3] = -3.0; 

      /* go */
      dummy = MPI_Wtime ();
      start = MPI_Wtime ();

      /* emmission d'un message au processus 1 */
      MPI_Isend (array, SIZE, MPI_DOUBLE,
		 dest, tag, MPI_COMM_WORLD,
		 &request);
      /* reception d'un message depuis le procesus 1 */
      MPI_Recv (array, SIZE, MPI_DOUBLE,
		dest, tagack, MPI_COMM_WORLD,
		&status);
      
      /* chrono */
      end = MPI_Wtime ();
      overhead = start - dummy ;
      elapsed = end - start - overhead ;

      /* memorisation des chronos */
      total += elapsed ;
      if (elapsed < min) min = elapsed ;
      if (elapsed > max) max = elapsed ; 
      
    } else if (self == dest) {
      /* ping ... pong */
      MPI_Recv (array, SIZE, MPI_DOUBLE,
		root, tag, MPI_COMM_WORLD,
		&status);
      array[SIZE-2] = -2.0 ; 
      MPI_Isend (array, SIZE, MPI_DOUBLE,
		 root, tagack, MPI_COMM_WORLD,
		 &request);
    }
  }

  /* verif integrite */
  if (self == root) {
    printf ("array[SIZE-2] = %g\n", array[SIZE-2]) ; 
    printf ("array[SIZE-3] = %g\n", array[SIZE-3]) ; 
  }
  
  /* resultat */
  if (self == root) {
    printf("%g Mo/s\n", (2*bytes*ITER)/total/1000000.);
    printf("(max %g Mo/s, min %g Mo/s)\n",
	   (2*bytes)/min/1000000.,
	   (2*bytes)/max/1000000.);
  }
  
  MPI_Finalize ();
}
