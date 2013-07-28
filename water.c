/**
 * Module: md.c
 * Author: Derek Williams
 * Created: 01/30/04
 * Purpose: 3-D hard spheres simulation of water molecules with periodic boundaries
 */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <GL/glut.h>

#include "list.h"
#include "water.h"

/**
 * Procedure: save
 * Purpose: save water molecule configuration as binary file
 * Inputs: name of file to be writted
 * Outputs: None
 **/
int save(char *fileName) 
{
     FILE *fp;
     fp = fopen(fileName, "w");
     if (!fp)
     {
      fprintf(stderr, "Error: opening %s for writing!\n", fileName);
      return 0;
     }

     fwrite(&props.density, sizeof(double), 1, fp);
	 fwrite(&props.rCut, sizeof(double), 1, fp);        // write out boundary dimensions
     fwrite(&props.deltaT, sizeof(double), 1, fp);       // write out step size
     fwrite(&stats.step, sizeof(int), 1, fp);            // write out step size
     fwrite(&numMols, sizeof(int), 1, fp);         // write number of molecules
     fwrite(mols, sizeof(waterMol), numMols, fp);  // write the molecules

     close(fp);
     return 1;
}

/**
 * Procedure: load
 * Purpose: load water molecule configuration from binary file
 * Inputs: name of file to be read
 * Outputs: None
 **/
int load(char *fileName)
{
     FILE *fp;
     int n;
     fp = fopen(fileName, "r");
     if (!fp) 
     {
      fprintf(stderr, "Error: opening %s for reading!\n", fileName);
      return 0;
     }

     fread(&props.density, sizeof(double), 1, fp);
     fread(&props.rCut, sizeof(double), 1, fp);        // write out boundary dimensions
     fread(&props.deltaT, sizeof(double), 1, fp); // read step size
     fread(&stats.step, sizeof(int), 1, fp);      // write out step size     
     fread(&numMols, sizeof(int), 1, fp);         // read number of molecules

     // now allocate space for molecules 
     mols = (waterMol*)malloc(sizeof(waterMol)*numMols);
//     printf("Allocating %d water molecules of size %d bytes\n", numMols, sizeof(waterMol));
//     printf("%d bytes total\n", sizeof(waterMol)*numMols);
     if (!mols) 
     {
      fprintf(stderr, "Unable to allocate molecules\n");
      return 0;
     }
     n = fread(mols, sizeof(waterMol), numMols, fp); // read the molecules
     if (n != numMols)
     {
      fprintf(stderr, "Error reading water molecules. Only %d of %d read\n", n, numMols);
      return 0;
     }
     distances = (double*)malloc(numMols*numMols*sizeof(double));
     if (!distances)
     {
      fprintf(stderr, "Unable to allocate distances\n");
      return 0;
     }

     close(fp);
     return 1;
}

/**
 * Procedure: openStatsFiles
 * Purpose: opens statistics files for writing statistics
 * Inputs: None
 * Outputs: None
 **/
void openStatsFiles(void)
{
     statfiles.energy = fopen(ENERGY_FILE, "w");
     statfiles.velocity = fopen(VELOCITY_FILE, "w");
     statfiles.momentum = fopen(MOMENTUM_FILE, "w");
     statfiles.angular = fopen(ANGULAR_FILE, "w");
     statfiles.pressuretemp = fopen(PRESSURETEMP_FILE, "w");
     statfiles.rdf = fopen(RDF_FILE, "w");
}

/**
 * Procedure: closeStatsFiles
 * Purpose: closes statistics files
 * Inputs: None
 * Outputs: None
 **/
void closeStatsFiles(void)
{
     close(statfiles.energy);
     close(statfiles.velocity);
     close(statfiles.momentum);
     close(statfiles.angular);
     close(statfiles.pressuretemp);
     close(statfiles.rdf);
}

/**
 * Procedure: doStats
 * Purpose: print out a line of statistical data to data files
 * Inputs: None
 * Outputs: None
 **/
void doStats(void)
{
//     printf("%d\n", stats.step);

     displayStats(stdout, statsType);

     displayStats(statfiles.energy, ENERGY_STATS);
     displayStats(statfiles.velocity, VELOCITY_X_STATS);
     displayStats(statfiles.momentum, MOMENTUM_STATS);
     displayStats(statfiles.angular, ANGULAR_VELOCITY_STATS);
     displayStats(statfiles.pressuretemp, PRESSURE_TEMPERATURE_STATS);
     displayStats(statfiles.rdf, RDF_STATS);
}

/**
 * Procedure: processArguments
 * Purpose: process the command line arguments
 * Inputs: the command line arguments
 * Outputs: success indication
 **/
int processArguments(int argc, char *argv[])
{
     statsType = 0;
     props.useGL = 0;

	 if (argc < 3)
     {
      printf("%s", helpStr);
      props.useGL = USE_GL_DEFAULT;
	  props.displayStep = DISPLAY_STEP_DEFAULT;
	  props.maxTime = MAX_TIME_DEFAULT;
     }
	 else
	 {
		 props.displayStep = atoi(argv[1]);
	     props.maxTime = atoi(argv[2]);
	 }
     if (argc >= 4)
      statsType = atoi(argv[3]);
     if (argc >= 5)
      props.useGL = atoi(argv[4]);
     if (argc == 6)
     {
      // initialize from previously saved data
      if (!load(argv[5])) return 0;
      initializePropsStats();
      props.equiliberationSteps = 0;
     }    
     else
     {
      // initialize velocities and positions
      if (!initialize()) 
	  {
		 printf("Unable to initialize!\n");
    	  return 0;
	  }
     }
     
     return 1;
}

/**
 * Main entry point of application
 * Inputs:
 *  argc - number of command line args
 *  argv - array of command line args
 */ 
int main(int argc, char *argv[])
{
     // initialize random number generator
     srand(time(NULL));
     siteMass[0] = 0.0;
     siteMass[1] = OXYGEN_MASS;
     siteMass[2] = HYDROGEN_MASS;
     siteMass[3] = HYDROGEN_MASS;

     openStatsFiles();

     if (!processArguments(argc, argv)) return 0;
     
     if (props.useGL)
     {
      // initialize the openGL runtime environment
      glutInit( &argc, argv );
      glutInitWindowSize(SIZE_X, SIZE_Y);
      glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
      glutCreateWindow(argv[0]);
      glutReshapeFunc(reshape);
      glutDisplayFunc(draw);
      glutIdleFunc(idle);

      initializeGL();

      // enter the infinite loop
      glutMainLoop();
     }
     else
     {
      // enter the infinite loop
      for(;;) idle();
     }

     // make the calling program happy :)
     return 0;
}

/**
 * Procedure: initializeGL
 * Purpose: initialize openGL environment
 * Inputs: None
 * Outputs: None
 **/
void initializeGL(void)
{
     GLfloat lin[] = { 0.7, 0.7, 0.7, 1.0 };
     GLfloat pos[] = { 0.0, 1.0, 1.0, 0.0 };
     double mh, mo;

     // initialize the clear color to white
     glClearColor(0.0, 0.0, 0.0, 0.0);

     glShadeModel(GL_SMOOTH);

     glEnable(GL_LIGHTING);
     glEnable(GL_LIGHT0);
     glEnable(GL_DEPTH_TEST);
     glEnable(GL_CULL_FACE);
     glEnable(GL_NORMALIZE);
     glEnable(GL_AUTO_NORMAL);
     
     glLightfv(GL_LIGHT0, GL_DIFFUSE, lin);
     glLightfv(GL_LIGHT0, GL_POSITION, pos);

     glFrontFace(GL_CCW);

     // create sphere list for rendereing of molecules
     props.oxygenSphere = glGenLists(1);
     glNewList(props.oxygenSphere, GL_COMPILE);
     glutSolidSphere(OXYGEN_RADIUS, 20, 20);
     glEndList();

     props.hydrogenSphere = glGenLists(1);
     glNewList(props.hydrogenSphere, GL_COMPILE);
     glutSolidSphere(HYDROGEN_RADIUS, 20, 20);
     glEndList();
}

/**
 * Procedure: 
 * Purpose: 
 * Inputs: None
 * Outputs: None
 **/
void assignCells(void)
{
     int i, x, y, z;
     nodeptr n;
     waterMol *m;

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];

      x = (int)floor(m->cm[0]/props.rCut);
      y = (int)floor(m->cm[1]/props.rCut);
      z = (int)floor(m->cm[2]/props.rCut);
     
      //assert(x < numCells[0] && y < numCells[1] && z < numCells[2]);

      m->cell[0] = x;
      m->cell[1] = y;
      m->cell[2] = z;

      n = createNode((void*)m);
      //assert(n);
      
      m->node = addNode(cells[x][y][z], n);
      //assert(m->node);
     }

//     printf("Cells assigned\n");
}

/**
 * Procedure: reassignCells
 * Purpose: reassign molecules to cells
 * Inputs: None
 * Outputs: None
 **/
void reassignCells(void)
{
     int i, x, y, z;
     nodeptr n;
     waterMol *m;

//     printf("Reassigning cells\n");

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];

      n = removeNode((nodeptr)m->node);
      //assert(n);

      x = (int)floor(m->cm[0]/props.rCut);
      y = (int)floor(m->cm[1]/props.rCut);
      z = (int)floor(m->cm[2]/props.rCut);

      //assert(x < numCells[0] && y < numCells[1] && z < numCells[2]);

//    printf("%d %d %d\n", x, y, z);

      m->cell[0] = x;
      m->cell[1] = y;
      m->cell[2] = z;
      
      m->node = addNode(cells[x][y][z], n);
      //assert(m->node);
     }

//     printf("Cells reassigned\n");
}

/**
 * Procedure: initializeCells
 * Purpose: initialize cells array for storing molecules
 * Inputs: None
 * Outputs: None
 **/
void initializeCells(void)
{
     int i, j, k;

     numCells[0] = (int)ceil(bounds[0][1]/props.rCut);
     numCells[1] = (int)ceil(bounds[1][1]/props.rCut);
     numCells[2] = (int)ceil(bounds[2][1]/props.rCut);

//     printf("cells (%d, %d, %d)\n", numCells[0], numCells[1], numCells[2]);

     // allocate the whole 3 dimensional cells array
     cells = (nodeptr***)malloc(numCells[0]*sizeof(nodeptr**));
     if (!cells)
     {
      fprintf(stderr, "Error: unable to allocate cells!\n");
      closeStatsFiles();
      exit(0);
     }
     for (i = 0; i < numCells[0]; i++)
     {
      cells[i] = (nodeptr**)malloc(numCells[1]*sizeof(nodeptr*));
      if (!cells) 
      {
           fprintf(stderr, "Error: unable to allocate cells!\n");
           closeStatsFiles();
           exit(0);
      }
      for (j = 0; j < numCells[1]; j++)
      {      
           cells[i][j] = (nodeptr*)malloc(numCells[2]*sizeof(nodeptr));
           if (!cells)
           {
            fprintf(stderr, "Error: unable to allocate cells!\n");
            closeStatsFiles();
            exit(0);
           }
           for (k = 0; k < numCells[2]; k++)
           {
            cells[i][j][k] = createNode(NULL);
            if (!cells[i][j][k])
            {
             fprintf(stderr, "Error: unable to create cell[%d][%d][%d]!\n", i, j, k);
             closeStatsFiles();
             exit(0);
            }
           }
      }
     }

//     printf("Cells initialized\n");
}

/**
 * Procedure: draw
 * Purpose: draw the profile display
 * Inputs: None
 * Outputs: None
 **/
static void draw(void)
{
     int i;   

     // clear the old drawing
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     glPushMatrix();

     glDisable(GL_LIGHTING);
     glDisable(GL_DEPTH_TEST);

//   draw the boundary
     glPushMatrix();
     glTranslated(0.0, 0.0, -bounds[2][1]);
     glutWireCube(bounds[0][1]);
     glPopMatrix();

     glTranslated(0.0, 0.0, -bounds[2][1]/2);

     glEnable(GL_LIGHTING);
     glEnable(GL_DEPTH_TEST);

     glTranslated(-bounds[0][1]/2, -bounds[1][1]/2, 0.0);

     for (i = 0; i < numMols; i++)
      drawMolecule(&mols[i]);

     glPopMatrix();

//   swap the buffers; force the polygon to be drawn on the screen
     glutSwapBuffers();
}

/**
 * Procedure: drawMolecule
 * Purpose: draw a 3D molecule onto the display
 * Inputs: the molecule to be drawn
 * Outputs: None
 **/
void drawMolecule(waterMol *mol)
{
     site *s1, *s2, *s3;
     GLfloat shn = 120;
     GLfloat oxygen[3][4] = { { 0.0, 0.0, 1.0, 0.0 },
                  { 0.01, 0.01, 1.0, 0.0 },
                  { 0.2, 0.2, 0.2, 0.0 } };
     GLfloat hydrogen[3][4] =   { { 0.1, 0.1, 0.1, 0.0 },
                  { 0.9, 0.9, 0.9, 0.0 },
                  { 0.2, 0.2, 0.2, 0.0 } };
      
     s1 = &mol->sites[1];
     s2 = &mol->sites[2];
     s3 = &mol->sites[3];

     glPushMatrix();
     glTranslated(mol->cm[0], mol->cm[1], -mol->cm[2]);
     
     // draw oxygen atom
     glMaterialfv(GL_FRONT, GL_AMBIENT, oxygen[0]);
     glMaterialfv(GL_FRONT, GL_DIFFUSE, oxygen[1]);
     glMaterialfv(GL_FRONT, GL_SPECULAR, oxygen[2]);
     glMaterialf(GL_FRONT, GL_SHININESS, shn);    
     
     glPushMatrix();
     glTranslated(s1->pos[0], s1->pos[1], -s1->pos[2]);
     glCallList(props.oxygenSphere);
     glPopMatrix();  

#if ROTATIONAL_MOTION
     // draw bonds
     glBegin(GL_LINES);
     glVertex3f(s1->pos[0], s1->pos[1], -s1->pos[2]);
     glVertex3f(s2->pos[0], s2->pos[1], -s2->pos[2]);  
     glVertex3f(s1->pos[0], s1->pos[1], -s1->pos[2]);
     glVertex3f(s3->pos[0], s3->pos[1], -s3->pos[2]);  
     glEnd();

     // draw hydrogen atoms
     glMaterialfv(GL_FRONT, GL_AMBIENT, hydrogen[0]);
     glMaterialfv(GL_FRONT, GL_DIFFUSE, hydrogen[1]);
     glMaterialfv(GL_FRONT, GL_SPECULAR, hydrogen[2]);
     glMaterialf(GL_FRONT, GL_SHININESS, shn);
     
     glPushMatrix();
     glTranslated(s2->pos[0], s2->pos[1], -s2->pos[2]);
     glCallList(props.hydrogenSphere);
     glPopMatrix();
     
     glPushMatrix();
     glTranslated(s3->pos[0], s3->pos[1], -s3->pos[2]);
     glCallList(props.hydrogenSphere);
     glPopMatrix();
#endif
     
     glPopMatrix();
}

/**
 * process window resize events
 * inputs:
 *  width, height - new width and height of window
 */ 
static void reshape(int w, int h)
{
     width = w;
     height = h;

     // setup our viewing area
     glViewport(0,0, (GLsizei) width, (GLsizei) height);

     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluPerspective(91, (double)width/(double)height, 1, 2*bounds[0][1]);
     glMatrixMode(GL_MODELVIEW);
     glLoadIdentity();
}

/**
 * Procedure: initializeVelocities
 * Purpose: initialize and scale the velocities of the molecues
 * Inputs: None
 * Outputs: None
 **/
void initializeVelocities(void)
{
     double avgvel[DIM];
     int i, j;
     waterMol *m;

     // initialize velocities random for now

//     printf("vMag = %lf\n", vMag);
     
     for (i = 0; i < DIM; i++)
      avgvel[i] = 0.0; 

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      randomUnitVector(m->veloc);
      for (j = 0; j < DIM; j++)
      {
           m->veloc[j] *= props.vMag;
           avgvel[j] += m->veloc[j];
      }
     }

     for (i = 0; i < DIM; i++)
      avgvel[i] /= (double)numMols; 

     // fix velocities so that total momentum is 0
     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      for (j = 0; j < DIM; j++)
           m->veloc[j] -= avgvel[j];
     }   
}

/**
 * Procedure: initializeAngularVelocities
 * Purpose: initialize and scale angular velocities of the molecules
 * Inputs: None
 * Outputs: None
 **/
void initializeAngularVelocities(void)
{
     double omega, e[3];
     int i, j;
     waterMol *m;

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      randomUnitVector(e);
      omega = 0.0;
      for (j = 0; j < DIM; j++)
           omega += m->inert[j] * SQR(e[j]);
      omega = props.vMag / sqrt(omega);
      m->qv[0] = 0.5 * omega * (m->q[3]*e[0] - m->q[2]*e[1] + m->q[1]*e[2]);
      m->qv[1] = 0.5 * omega * (m->q[2]*e[0] + m->q[3]*e[1] - m->q[0]*e[2]);
      m->qv[2] = 0.5 * omega * (m->q[1]*e[0] + m->q[0]*e[1] + m->q[3]*e[2]);
      m->qv[3] = 0.5 * omega * (m->q[0]*e[0] - m->q[1]*e[1] - m->q[2]*e[2]);
     }
}

/**
 * Procedure: predict
 * Purpose: 
 * Inputs: None
 * Outputs: None
 **/
void predict(void)
{
#if USE_GEAR
     predictGear();
#else
     predictVerlet();
#endif     
}

/**
 * Procedure: correct
 * Purpose: 
 * Inputs: None
 * Outputs: None
 **/
void correct(void)
{
#if USE_GEAR
     correctGear();
#else
     correctVerlet();
#endif     
}

/**
 * Procedure: predictRotations
 * Purpose: 
 * Inputs: None
 * Outputs: None
 **/
void predictRotations(void)
{
#if USE_GEAR
     predictRotationsGear();
#else
     predictRotationsVerlet();
#endif     
}

/**
 * Procedure: correctRotations
 * Purpose: 
 * Inputs: None
 * Outputs: None
 **/
void correctRotations(void)
{
#if USE_GEAR
     correctRotationsGear();
#else
     correctRotationsVerlet();
#endif     
}

/**
 * Procedure: idle
 * Purpose: do deltaT time steps of the algorithm
 * Inputs: None
 * Outputs: None
 **/
static void idle(void)
{
     int i, j, x, y;
     time_t currentTime;

     predict();
     predictRotations();
#if !ALL_PAIRS
     reassignCells();
#endif
     evaluateForces();
     computeRotationalAccelerations();
     correct();
     correctRotations();

     adjustTemperature();
     updateStats();

     if ((stats.step % props.displayStep) == 0)
     {
      doStats();      
      if (props.useGL) glutPostRedisplay();
     }
     if (props.equiliberationSteps)
     {
      if (stats.step == props.equiliberationSteps)
           resetAccumulators();
     }
     if (stats.step >= props.maxTime)
     {
      save(DATA_FILE_NAME);  // save the data	  
	  closeStatsFiles();
	  exit(EXIT_SUCCESS);
     }

     stats.time += props.deltaT;
     stats.step ++;
}

/**
 * Procedure: resetAccumulators
 * Purpose: reset accumulators used for computing time canonical statistics
 * Inputs: None
 * Outputs: None
 **/
void resetAccumulators(void)
{
     int i;

     stats.accN = 1;
     stats.nRDF = 0;
     stats.accKE = 0.0;     
     stats.accAKE = 0.0;     
     stats.accVIRIAL = 0.0;

     for (i = 0; i < RDF_BINS; i++)
     {
      stats.accRDFOO[i] = 0;
      stats.accRDFHH[i] = 0;
      stats.accRDFOH[i] = 0;
     }
}

/**
 * Procedure: updateStats
 * Purpose: update the statistical data
 * Inputs: None
 * Outputs: None
 **/
void updateStats(void)
{
     stats.accKE += stats.KE; 
     stats.accAKE += stats.AKE;     
     stats.accVIRIAL += stats.VIRIAL;
     stats.meanKE = stats.accKE/stats.accN; 
     stats.meanAKE = stats.accAKE/stats.accN;
#if ROTATIONAL_MOTION
     stats.T = 2.0/(double)(DIM*numMols)*(stats.meanKE + stats.meanAKE);
#else
     stats.T = 2.0/(double)(DIM*numMols)*stats.meanKE;
#endif
/*
#if ROTATIONAL_MOTION
     stats.T = 2.0*(stats.KE + stats.AKE)/(double)(DIM*numMols);
#else
     stats.T = 2.0*stats.KE/(double)(DIM*numMols);
#endif
*/
     stats.P = numMols*stats.T + stats.accVIRIAL/(double)(DIM*stats.accN); 
     stats.P /= stats.V;
     stats.accN ++;
}

/**
 * Procedure: initialize
 * Purpose: initialize the algorithm
 * Inputs: None
 * Outputs: None
 **/
int initialize(void)
{
     int i, j;     
     waterMol *m;
   
     numMols = NUM_MOLECULES;
     mols = (waterMol*)malloc(sizeof(waterMol)*numMols);
     if (!mols)
     {
      fprintf(stderr, "Unable to allocate molecules\n");
      return 0;
     }
     distances = (double*)malloc(numMols*numMols*sizeof(double));
     if (!distances)
     {
      fprintf(stderr, "Unable to allocate distances\n");
      return 0;
     }

     for (i = 0; i < numMols; i++)
     {
      initializeMolecule(&mols[i]);
      mols[i].n = i;
     }

     props.density = DENSITY;
     props.rCut = CUTOFF_DISTANCE;
     props.equiliberationSteps = EQUILIBERATION_TIME;

     initializePropsStats();

     buildLattice();
     initializeVelocities();
#if ROTATIONAL_MOTION
     initializeAngularVelocities();
#endif
#if !ALL_PAIRS
     initializeCells();
     assignCells();
#endif
	return 1;
}

/**
 * Procedure: computeKE
 * Purpose: calculate the total kinetic energy of the isolated system
 * Inputs: None
 * Outputs: None
 **/
void computeKE(void)
{
     int i, j;
     waterMol *m;

     // computer kinetic energy and linear momentum of the system assuming masses are 1 for now
     stats.KE = 0.0;
     stats.AKE = 0.0;
     stats.LM[0] = stats.LM[1] = stats.LM[2] = 0.0;
     stats.AM[0] = stats.AM[1] = stats.AM[2] = 0.0;
     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      for (j = 0; j < DIM; j++)
      {
           stats.KE += m->veloc[j] * m->veloc[j];
           stats.LM[j] += m->veloc[j];
           stats.AKE += m->omega[j] * m->omega[j];
           stats.AM[j] += m->omega[j];
      }
     }
     stats.KE *= 0.5; 
     stats.AKE *= 0.5; 
}

/**
 * Procedure: initializePropsStats
 * Purpose: initialize general properties of the system
 * Inputs: None
 * Outputs: None
 **/
void initializePropsStats(void)
{
     int I;     
     double l, a;

     l = cbrt((double)numMols/props.density);

     bounds[0][0] = bounds[1][0] = bounds[2][0] = 0.0;
     bounds[0][1] = bounds[1][1] = bounds[2][1] = l;

     stats.V = l*l*l;

     stats.rdfMax = 0.5*l;
     stats.dr = stats.rdfMax/(double)RDF_BINS;

     props.vMag = sqrt((double)DIM * (1.0 - 1.0 / (double)(numMols)) * TEMPERATURE);
#if ROTATIONAL_MOTION
     props.vMag /= 2.0; // half goes to translation and half to kinetic
#endif
#if FIXED_SITES
     props.vMag = 0.0;
#endif

     stats.time = 0;
     stats.step = 1;

     resetAccumulators();

     time(&start);

     props.deltaT = DT;
     props.deltaT2 = DT*DT;
     props.tempAdjustStep = ADJUST_TEMP_STEP;

     props.rrCut = props.rCut*props.rCut;

//     props.sigma = pow(LJ_A/LJ_C, 1.0/6.0); 
//     props.epsilon = LJ_C/(4.0*pow(props.sigma, 6.0));
}

/**
 * Procedure: buildLattice
 * Purpose: build cubic lattice on which to put the molecules
 * Inputs: None
 * Outputs: None
 **/
void buildLattice(void)
{
     double a, r[4][DIM], I;
     int i, m, x, y, z;     

     I = cbrt((double)numMols/4.0);
     a = bounds[0][1] / I;

     // initialize lattice site vectors
     r[0][0] = 0.0;    r[0][1] = 0.0;   r[0][2] = 0.0; 
     r[1][0] = 0.0;    r[1][1] = a/2.0; r[1][2] = a/2.0; 
     r[2][0] = a/2.0;  r[2][1] = 0.0;   r[2][2] = a/2.0; 
     r[3][0] = a/2.0;  r[3][1] = a/2.0; r[3][2] = 0.0; 

     m = 0;
     for (x = 0; x < I; x++)
     {
      for (y = 0; y < I; y++)
      {
           for (z = 0; z < I; z++)
           {      
            for (i = 0; i < 4; i++, m++)
            {
             mols[m].cm[0] = x*a + r[i][0]; 
             mols[m].cm[1] = y*a + r[i][1]; 
             mols[m].cm[2] = z*a + r[i][2];
             //assert(m < numMols);
             fixPosition(m);
            }
           }
      }
     }
     //assert(m == numMols);
}

/**
 * Procedure: randomUnitVector
 * Purpose: make a normalized random vector
 * Inputs: the vector
 * Outputs: None
 **/
void randomUnitVector(double vec[DIM])
{
     int i;
     double mag, sign, rnd;
     mag = 0.0;
     for (i = 0; i < DIM; i++)
     {
      sign = (rand()/(double)RAND_MAX < 0.5) ? -1.0 : 1.0; 
      vec[i] = sign * rand()/(double)RAND_MAX;
      mag += vec[i]*vec[i];
     }
     mag = sqrt(mag);
     for (i = 0; i < DIM; i++)
      vec[i] /= mag;
}

/**
 * Procedure: randomQuarternion
 * Purpose: make a normalized random quarternion coordinate
 * Inputs: the quarternion
 * Outputs: None
 **/
void randomQuarternion(double q[4])
{
     int i;
     double sign, mag;
     mag = 0.0;
     for (i = 0; i < 4; i ++)
     {
      sign = (rand()/(double)RAND_MAX < 0.5) ? 1.0 : -1.0;
      q[i] = sign * rand() / (double)RAND_MAX;
     }
     normalizeQuarternion(q);
}

/**
 * Procedure: normalizeQuarternion
 * Purpose: normalize a quarternion
 * Inputs: the quarternion
 * Outputs: None
 **/
void normalizeQuarternion(double q[4])
{
     int i;
     double mag;
     mag = 0.0;
     for (i = 0; i < 4; i ++)
      mag += q[i]*q[i];
     for (i = 0; i < 4; i++)
      q[i] /= sqrt(mag); 
}

/**
 * Procedure: computeInertia
 * Purpose: comput the moment of inertia of the molecule
 * Inputs: the molecule
 * Outputs: None
 **/
void computeInertia(waterMol *mol)
{
     int i;
     site *s;

     // compute moment of inertia
     for (i = 0; i < NUM_SITES; i++)
     {
      s = &mol->sites[i];

      mol->inert[0] += siteMass[i] * (s->pos[1] * s->pos[1] + s->pos[2] * s->pos[2]);
      mol->inert[1] += siteMass[i] * (s->pos[0] * s->pos[0] + s->pos[2] * s->pos[2]);
      mol->inert[2] += siteMass[i] * (s->pos[0] * s->pos[0] + s->pos[1] * s->pos[1]);
     }
}

/**
 * Procedure: initializeMolecule
 * Purpose: initialize/create a molecule (molecule must be allocated!)
 * Inputs: the molecule
 * Outputs: None
 **/
void initializeMolecule(waterMol *mol)
{
     int i;
     double len, pos[DIM];
     site *s;
     for (i = 0; i < DIM; i++)
     {
      mol->cm[i] = 0.0;
      mol->veloc[i] = 0.0;
      mol->cmold[i] = 0.0;
      mol->velocold[i] = 0.0;
      mol->accel[i] = 0.0;
      mol->accel1[i] = 0.0;
      mol->accel2[i] = 0.0;
      mol->euler[i] = 0.0;
      mol->omega[i] = 0.0;
      mol->tau[i] = 0.0;
      mol->alpha[i] = 0.0;
     }
     setSites(mol);
//     computeInertia(mol);
     mol->inert[0] = 0.0098;
     mol->inert[1] = 0.0034;
     mol->inert[2] = 0.0064;

     for (i = 0; i < 4; i++)
     {
      mol->qv[i] = mol->qa[i] = 0.0;
      mol->qvold[i] = mol->qold[i] = 0.0;
      mol->qa1[i] = mol->qa2[i] = 0.0;
     }

     randomQuarternion(mol->q);

//     rotateMolecule(mol, euler);
     for (i = 0; i < NUM_SITES; i++)
      rotateToBodyFrame(mol->sites[i].pos, mol->q);
}

/**
 * Procedure: setSites
 * Purpose: setup the atomic interaction sites of a molecule
 * Inputs: a molecule
 * Outputs: None
 **/
void setSites(waterMol *mol)
{
//     double xh, yh, ycm;

//     xh = R_OH * COS( 0.5 * THETA_HOH );
//     yh = R_OH * SIN( 0.5 * THETA_HOH );

//     ycm = (xh + yh) / 18.0;

     site *s;
     s = &mol->sites[0];
     s->type = M_SITE;
     s->pos[0] = 0.0;
     s->pos[1] = 0.0;
     s->pos[2] = 0.0274;

     s = &mol->sites[1];
     s->type = O_SITE;
     s->pos[0] = 0.0;
     s->pos[1] = 0.0;
     s->pos[2] = -0.0206;

     s = &mol->sites[2];
     s->type = H_SITE;
     s->pos[0] = 0.0;
     s->pos[1] = 0.240;
     s->pos[2] = 0.165;

     s = &mol->sites[3];
     s->type = H_SITE;
     s->pos[0] = 0.0;
     s->pos[1] = -0.240;
     s->pos[2] = 0.165;
}

/**
 * Procedure: displayStats
 * Purpose: format statistics for display and display them
 * Inputs: file to display to and bit flag indicating type(s) of statistics to display
 * Outputs: None
 **/
void displayStats(FILE* fp, int type)
{
     int i, bin;
     double vel, P, T, meanKE, norm, V, goo, goh, ghh, rinner, router, rho;
     time_t now;
     switch (type)
     {
     case ENERGY_STATS:
#if ROTATIONAL_MOTION
      fprintf(fp, "%d\t%.4lf\t%.4lf\n", stats.step, stats.KE/(double)numMols + stats.AKE/(double)numMols, 
          stats.UE/(double)numMols);  
#else
      fprintf(fp, "%d\t%.4lf\t%.4lf\n", stats.step, stats.KE/(double)numMols, stats.UE/(double)numMols);  
#endif
      break;
     case MOMENTUM_STATS:
      fprintf(fp, "%d\t%.2le\t%.2le\t%.2le\t%.2le\t%.2le\t%.2le\n", stats.step, stats.LM[0], stats.LM[1], 
          stats.LM[2], stats.AM[0], stats.AM[1], stats.AM[2]);  
      break;
     case VELOCITY_X_STATS:
      for (i = 0; i < numMols; i++)
           fprintf(fp,"%.2le ", mols[i].veloc[0]);
      fprintf(fp, "\n");
      break;
     case VELOCITY_Y_STATS:
      for (i = 0; i < numMols; i++)
           fprintf(fp, "%.2le ", mols[i].veloc[1]);
      fprintf(fp, "\n");
      break;
     case VELOCITY_Z_STATS:
      for (i = 0; i < numMols; i++)
           fprintf(fp, "%.2le ", mols[i].veloc[2]);
      fprintf(fp, "\n");
      break;
     case ANGULAR_VELOCITY_STATS:
      for (i = 0; i < numMols; i++)
      {
           vel = sqrt(mols[i].omega[0] * mols[i].omega[0] + 
              mols[i].omega[1] * mols[i].omega[1] + 
              mols[i].omega[2] * mols[i].omega[2]);
           fprintf(fp, "%.2le ", vel);
      }
      fprintf(fp, "\n");
      break;
     case FPS_STATS:
      time(&now);
      fprintf(fp, "%lf frames/sec\n", (double)stats.step/((double)now-(double)start));
      break;
     case MOLECULE_STATS:
      for (i = 0; i < numMols; i++)
      {
           fprintf(fp, "Molecule %d:\n", i);
           fprintf(fp, "\tcm (%lf, %lf, %lf)\n", mols[i].cm[0], mols[i].cm[1], mols[i].cm[2]);
           fprintf(fp, "\tveloc (%lf, %lf, %lf)\n", mols[i].veloc[0], mols[i].veloc[1], mols[i].veloc[2]);
           fprintf(fp, "\taccel (%lf, %lf, %lf)\n", mols[i].accel[0], mols[i].accel[1], mols[i].accel[2]);
           fprintf(fp, "\n");
      }
      break;
     case PRESSURE_TEMPERATURE_STATS:
      fprintf(fp, "%d\t%.2lf\t%.2lf\n", stats.step, stats.P, stats.T);
      break;
     case QUATERNION_STATS:
      for (i = 0; i < numMols; i++)
           fprintf(fp, "(%.2lf, %.2lf, %.2lf, %.2lf) ", mols[i].q[0], mols[i].q[1], mols[i].q[2], mols[i].q[3]);
      fprintf(fp, "\n");
      break;    
     case PAIR_DISTRIBUTION_STATS:
      for (i = 0; i < pairs; i++)
           fprintf(fp, "%lf ", distances[i]);
      fprintf(fp, "\n");
      break;
     case RDF_STATS:
      if ((stats.step % props.maxTime) == 0)
      {
           rho = numMols/stats.V;
//         printf("rho = %lf\tnRDF = %lf\n", rho, stats.nRDF);
           norm = (double)numMols*stats.nRDF*rho;

           for (bin = 0; bin < RDF_BINS; bin++)
           {
            rinner = (double)bin*stats.dr;
            router = rinner + stats.dr;

            V = 4.0/3.0*M_PI*(pow(router, 3.0) - pow(rinner, 3.0));
//          printf("%d\t%lf\t%lf\t%lf\n", bin, rinner, router, V);
            goo = (double)(stats.accRDFOO[bin])/norm/V;
            ghh = (double)(stats.accRDFHH[bin])/norm/V/4.0;
            goh = (double)(stats.accRDFOH[bin])/norm/V/2.0;
            fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", rinner+0.5*stats.dr, goo, ghh, goh);
           }
      }
      break;
     }
}

/**
 * Procedure: adjustTemperature
 * Purpose: adjust the temperature of the system so that it reaches desired temperature
 * Inputs: None
 * Outputs: None
 **/
void adjustTemperature(void)
{
     int i, j, bin;
     double fac, ke, ake;
     waterMol *m;

     if (props.tempAdjustStep)// && stats.step < props.equiliberationSteps)
     {
      if ((stats.step % props.tempAdjustStep) == 0)
      {
//         ke = stats.meanKE/(double)numMols;
           ke = stats.KE/(double)numMols;
           fac = props.vMag/sqrt(2.0*ke);
           for (i = 0; i < numMols; i++)
           {
            m = &mols[i];
            for (j = 0; j < DIM; j++)
             m->veloc[j] *= fac;
           }
#if ROTATIONAL_MOTION
//               ake = stats.meanAKE/(double)numMols;
           ake = stats.AKE/(double)numMols;
           fac = props.vMag/sqrt(2.0*ake);
           for (i = 0; i < numMols; i++)
           {
            m = &mols[i];
            for (j = 0; j < DIM; j++)
             m->qv[j] *= fac;
           }
#endif 
      }
     }
}

/**
 * Procedure: evaluateForces
 * Purpose: calculate the forces acting at each interaction site
 * Inputs: None
 * Outputs: None
 **/
void evaluateForces(void)
{
     int i, j, k, x, y, z, dx, dy, dz, nx, ny, nz, xx, yy, zz, n, bin;
     double u, F[DIM], r2, dr[DIM], v, r;
     nodeptr *next, *tmp;
     waterMol *m1, *m2;
     char s[100];

     stats.UE = u = 0.0;
     stats.VIRIAL = v = 0.0;

     // initialize accelerations and site forces
     for (i = 0; i < numMols; i++)
     {
      for (j = 0; j < DIM; j++)
      {
           mols[i].accel[j] = 0.0;
           mols[i].alpha[j] = 0.0;
      }

      setSites(&mols[i]);

      // initialize site forces to 0
      for (j = 0; j < NUM_SITES; j++)
      {
           for (k = 0; k < DIM; k++)
            mols[i].sites[j].force[k] = 0.0;
           rotateToSpaceFrame(mols[i].sites[j].pos, mols[i].q);
      }   
     }     

     n=0;

#if ALL_PAIRS

     // compute the total forces on each of the particles
     for (i = 0, k = 0; i < numMols - 1; i++)
     {
      for (j = i + 1; j < numMols; j++, k++)
      {
           // compute the center to center distance between the two molecules
           r2 = distance(&mols[i], &mols[j], dr);
           
           if ((stats.step % props.displayStep) == 0)
           {
            // o-o rdf
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[1], &mols[j].sites[1], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFOO[bin] += 1;
            }
            // o-h rdf
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[1], &mols[j].sites[2], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFOH[bin] += 1;
            }
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[1], &mols[j].sites[3], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFOH[bin] += 1;
            }
            // h-h rdf
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[2], &mols[j].sites[2], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFHH[bin] += 1;
            }
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[3], &mols[j].sites[3], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFHH[bin] += 1;
            }
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[3], &mols[j].sites[2], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFHH[bin] += 1;
            }
            r = sqrt(siteDistance(&mols[i], &mols[j], &mols[i].sites[2], &mols[j].sites[3], dr));
                if (r < stats.rdfMax)
            {
             bin = (int)floor( r / stats.dr );
             stats.accRDFHH[bin] += 1;
            }
           }

           if (r2 < props.rrCut) 
           {
//          printf("r2 = %lf\n", r2);
            computeSiteForces(&mols[i], &mols[j], &u, &v);
            stats.UE += u; // accumulate potential energy
            stats.VIRIAL += v;
            n++;
           }
      }
     }
     if ((stats.step % props.displayStep) == 0)
      stats.nRDF ++;
     pairs = k;
#else
     nx = numCells[0]; ny = numCells[1]; nz = numCells[2];
     for (i = 0, k = 0; i < numMols; i++)
     {
      m1 = &mols[i];
      x = m1->cell[0]; y = m1->cell[1]; z = m1->cell[2];
//    printf("Scanning molecule %d of %d in cell (%d, %d, %d)...\n", i, numMols, x, y, z);

//    printNode(m1->node, "");

      tmp = removeNode(m1->node);
      //assert(tmp); // don't calculate force on itself
//    sprintf(s, "%d", i);
//    printNode(tmp, s);

      for (dx = -1; dx <= 1; dx++)
      {
           for (dy = -1; dy <= 1; dy++)
           {
            for (dz = -1; dz <= 1; dz++)
            {
             xx = x+dx; yy = y+dy; zz = z+dz;
//           printf("(%d, %d, %d) ", xx, yy, zz);
             if (xx >= nx) xx = 0; if (xx < 0) xx = nx-1;
             if (yy >= ny) yy = 0; if (yy < 0) yy = ny-1;
             if (zz >= nz) zz = 0; if (zz < 0) zz = nz-1;
//           printf("(%d, %d, %d) ", xx, yy, zz);
             //assert(xx < numCells[0] && yy < numCells[1] && zz < numCells[2]);
             //assert(cells[xx][yy][zz]);
             next = cells[xx][yy][zz]->next;
             while (next)
             {
                  m2 = (waterMol*)next->data;
                  //assert(m2->n >= m1->n);
                  // compute the center to center distance between the two molecules
                  r2 = distance(m1, m2, dr);
                  distances[k] = r2;
                  k++;

                  if (r2 < props.rrCut) 
                  {
                   computeSiteForces(m1, m2, &u, &v);
                   stats.UE += u; // accumulate potential energy
                   stats.VIRIAL += v;
                   n++;
                  }

                  next = next->next;
             }
            }
           }
      }
//    tmp = addNode(cells[x][y][z], tmp);
//    //assert(tmp);
//    printf("\nScanning complete\n", i, x, y, z);
     }
     pairs = k;
#endif
//     printf("n = %d\n", n);
}

/**
 * Procedure: computeRotationalAccelerations
 * Purpose: compute rotational acceleration of each molecule
 * Inputs: None
 * Outputs: None
 **/
void computeRotationalAccelerations(void)
{
#if ROTATIONAL_MOTION
     int i, j, k;
     waterMol *m;
     double a[4], WT[4][4];

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];

      computeAngularVelocity(m);
      computeTorque(m);

      a[0] = (m->tau[0] + (m->inert[1] - m->inert[2]) * m->omega[1] *  m->omega[2]) / m->inert[0];
      a[1] = (m->tau[1] + (m->inert[2] - m->inert[0]) * m->omega[0] *  m->omega[2]) / m->inert[1];
      a[2] = (m->tau[2] + (m->inert[0] - m->inert[1]) * m->omega[0] *  m->omega[1]) / m->inert[2];
      a[3] = 0.0; for (j = 0; j < 4; j++) a[3] -= 2.0 * m->qv[j]*m->qv[j];

      WT[0][0] = W11(m->q);  WT[0][1] = W21(m->q);  WT[0][2] = W31(m->q);  WT[0][3] = W41(m->q);
      WT[1][0] = W12(m->q);  WT[1][1] = W22(m->q);  WT[1][2] = W32(m->q);  WT[1][3] = W42(m->q);
      WT[2][0] = W13(m->q);  WT[2][1] = W23(m->q);  WT[2][2] = W33(m->q);  WT[2][3] = W43(m->q);
      WT[3][0] = W14(m->q);  WT[3][1] = W24(m->q);  WT[3][2] = W34(m->q);  WT[3][3] = W44(m->q);

      for (j = 0; j < 4; j++) 
           m->qa[j] = 0.0;
      for (j = 0; j < 4; j++) 
      {
           for (k = 0; k < 4; k++)
            m->qa[j] += 0.5 * WT[j][k] * a[k];
      } 
     }
#endif
}

/**
 * Procedure: rotateToBodyFrame
 * Purpose: rotate vector from space to body frame of reference
 * Inputs: the vector and the quarternion
 * Outputs: None
 **/
void rotateToBodyFrame(double v[DIM], double q[4])
{
     int i, j;
     double A[3][3];
     double x[3];

     // copy the vector
     for (i = 0; i < DIM; i++)
     {
      x[i] = v[i];
      v[i] = 0.0;
     }

     // create rotation matrix 
     A[0][0] = A11(q); A[0][1] = A12(q); A[0][2] = A13(q);
     A[1][0] = A21(q); A[1][1] = A22(q); A[1][2] = A23(q);
     A[2][0] = A31(q); A[2][1] = A32(q); A[2][2] = A33(q);

     // multiply the vector by the matrix
     for (i = 0; i < DIM; i++)  
     {
      for (j = 0; j < DIM; j++)
           v[i] += 2.0*A[i][j]*x[j];
     }
}

/**
 * Procedure: 
 * Purpose: rotate vector from body to space frame of reference
 * Inputs: the vector and the quarternion
 * Outputs: None
 **/
void rotateToSpaceFrame(double v[DIM], double q[4])
{
     int i, j;
     double A[3][3];
     double x[3];

     // copy the vector
     for (i = 0; i < DIM; i++)
     {
      x[i] = v[i];
      v[i] = 0.0;
     }

     // create transpose rotation matrix 
     A[0][0] = A11(q); A[0][1] = A21(q); A[0][2] = A31(q);
     A[1][0] = A12(q); A[1][1] = A22(q); A[1][2] = A32(q);
     A[2][0] = A13(q); A[2][1] = A23(q); A[2][2] = A33(q);

     // multiply the vector by the matrix
     for (i = 0; i < DIM; i++)  
     {
      for (j = 0; j < DIM; j++)
           v[i] += 2.0*A[i][j]*x[j];
     }
}

/**
 * Procedure: computeTorque
 * Purpose: calculate the torque acting on a molecule
 * Inputs: the molecule
 * Outputs: None
 **/
void computeTorque(waterMol *m)
{
     double inert[DIM];
     int i, j, k;
     site *s;

     for (i = 0; i < DIM; i++)
      m->tau[i] = 0.0;

     for (i = 0; i < NUM_SITES; i++)
     {
      s = &m->sites[i];

      // use cross product to compute torque
      m->tau[0] += s->pos[1] * s->force[2] - s->pos[2] * s->force[1];
      m->tau[1] += s->pos[2] * s->force[0] - s->pos[0] * s->force[2];
      m->tau[2] += s->pos[0] * s->force[1] - s->pos[1] * s->force[0];

//    printf("%.2lf, %.2lf, %.2lf\n", m->tau[0], m->tau[1], m->tau[2]);
     }

//     printf("->%.2lf, %.2lf, %.2lf\n", m->tau[0], m->tau[1], m->tau[2]);

     rotateToBodyFrame(m->tau, m->q);
}

/**
 * Procedure: computeAngularVelocity
 * Purpose: calculate the angular velocity of a molecule
 * Inputs: the molecule
 * Outputs: None
 **/
void computeAngularVelocity(waterMol *m)
{
     int i, j;
     double W[4][4];

#define Q m->q
     W[0][0] = Q[3];  W[0][1] = Q[2];  W[0][2] =-Q[1];  W[0][3] =-Q[0];
     W[1][0] =-Q[2];  W[1][1] = Q[3];  W[1][2] = Q[0];  W[1][3] =-Q[1];
     W[2][0] = Q[1];  W[2][1] =-Q[0];  W[2][2] = Q[3];  W[2][3] =-Q[2];
     W[3][0] = Q[0];  W[3][1] = Q[1];  W[3][2] = Q[2];  W[3][3] = Q[3];
#undef Q

     for (i = 0; i < 3; i++) 
      m->omega[i] = 0.0;

     for (i = 0; i < 3; i++) 
     {
      for (j = 0; j < 4; j++)
           m->omega[i] += 2.0 * W[i][j] * m->qv[j];
     }
}

/**
 * Procedure: fixPosition
 * Purpose: apply perodic boundary conditions on x, y, and z to a molecule 
 *			constraining its position to the boundaries of space
 * Inputs: array index of molecule
 * Outputs: None
 **/
void fixPosition(int i)
{
     int j;
     for (j = 0; j < DIM; j++)
     {
      // periodic on x
      if (mols[i].cm[j] < bounds[j][0])
           mols[i].cm[j] += bounds[j][1];
      else if (mols[i].cm[j] >= bounds[j][1])
           mols[i].cm[j] -= bounds[j][1];
     }
}

/**
 * Procedure: rotateEuler
 * Purpose: perform euler rotation on a vector
 * Inputs: the euler angle and the vector
 * Outputs: None
 **/
void rotateEuler(double euler[DIM], double v[DIM])
{
     int i, j;
     double A[3][3];
     double x[3], cphi, ctheta, cpsi, sphi, stheta, spsi;

     // copy the vector
     for (i = 0; i < DIM; i++)
     {
      x[i] = v[i];
      v[i] = 0.0;
     }

     // calculate the sines and cosines
     cphi = cos(euler[0]); ctheta = cos(euler[1]); cpsi = cos(euler[2]);
     sphi = sin(euler[0]); stheta = sin(euler[1]); spsi = sin(euler[2]);

     // create inverse rotation matrix (pitch, roll, yaw convention)
     A[0][0] = ctheta*cphi; A[0][1] = spsi*stheta*cphi - cpsi*sphi; A[0][2] = cpsi*stheta*cphi + spsi*sphi;
     A[1][0] = ctheta*sphi; A[1][1] = spsi*stheta*sphi + cpsi*cphi; A[1][2] = cpsi*stheta*sphi - spsi*cphi;
     A[2][0] = -stheta;     A[2][1] = ctheta*spsi;                  A[2][2] = ctheta*cpsi;

     // multiply the vector by the matrix
     for (i = 0; i < DIM; i++)  
     {
      for (j = 0; j < DIM; j++)
           v[i] += A[i][j]*x[j];
     }
}

/**
 * Procedure: predictVerlet
 * Purpose: perform verlet predictor step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void predictVerlet(void)
{
     int i, j;
     waterMol *m;

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      // update positions
      for (j = 0; j < DIM; j++)
      {
           m->cm[j] += props.deltaT*m->veloc[j];
           m->cm[j] += 0.5*props.deltaT2*m->accel[j]; 
           m->veloc[j] += 0.5*props.deltaT*m->accel[j];

           //assert(m->cm[j] < 2*bounds[j][1]);
           //assert(m->cm[j] > -2*bounds[j][1]);      
           //assert(m->veloc[j] < MAX_VELOC);
      }

      fixPosition(i);
     }
}

/**
 * Procedure: correctVerlet
 * Purpose: perform verlet corrector step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void correctVerlet(void)
{
     int i, j;
     waterMol *m;

     stats.KE = 0;
     stats.LM[0] = stats.LM[1] = stats.LM[2] = 0;
     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      // update velocities
      for (j = 0; j < DIM; j++)
      {
           m->veloc[j] += 0.5*props.deltaT*m->accel[j]; // correct velocity     

           stats.KE += m->veloc[j]*m->veloc[j];                      // kinetic energy
           stats.LM[j] += m->veloc[j];                                    // linear momentum
      }   
      fixPosition(i);
     }
     stats.KE *= 0.5;
}

/**
 * Procedure: predictGear
 * Purpose: perform gear predictor step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void predictGear(void)
{
     int i, j;
     waterMol *m;
     double cr[] = {19.0, -10.0, 3.0};
     double cv[] = {27.0, -22.0, 7.0};
     double div = 24.0;

     stats.KE = 0;
     stats.LM[0] = stats.LM[1] = stats.LM[2] = 0;
     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      for (j = 0; j < DIM; j++)
      {
           m->cmold[j] = m->cm[j];
           m->velocold[j] = m->veloc[j];

               // update positions, velocities, accelerations, and other time derivatives
           m->cm[j] += props.deltaT*m->veloc[j];
           m->cm[j] += (props.deltaT2/div) * (cr[0]*m->accel[j] + cr[1]*m->accel1[j] + cr[2]*m->accel2[j]); 
           m->veloc[j] = (m->cm[j] - m->cmold[j])/props.deltaT;
           m->veloc[j] += (props.deltaT/div) * (cv[0]*m->accel[j] + cv[1]*m->accel1[j] + cv[2]*m->accel2[j]);
           m->accel2[j] = m->accel1[j];
           m->accel1[j] = m->accel[j];

           // update kinetic energy and momentum
           stats.KE += m->veloc[j]*m->veloc[j];
           stats.LM[j] += m->veloc[j];         
      }
      fixPosition(i);
     }
     stats.KE *= 0.5;
}

/**
 * Procedure: correctGear
 * Purpose: perform gear corrector step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void correctGear(void)
{
     int i, j;
     waterMol *m;
     double cr[] = {3.0, 10.0, -1.0};
     double cv[] = {7.0, 6.0, -1.0};
     double div = 24.0;

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];
      for (j = 0; j < DIM; j++)
      {
           // correct positions and velocities
           m->cm[j] = m->cmold[j] + props.deltaT*m->velocold[j];
           m->cm[j] += (props.deltaT2/div) * (cr[0]*m->accel[j] + cr[1]*m->accel1[j] + cr[2]*m->accel2[j]); 
           m->veloc[j] = (m->cm[j] - m->cmold[j])/props.deltaT;
           m->veloc[j] += (props.deltaT/div) * (cv[0]*m->accel[j] + cv[1]*m->accel1[j] + cv[2]*m->accel2[j]);
      
//         //assert(mols[i].cm[j] < 2*bounds[j][1]);
//         //assert(mols[i].cm[j] > -bounds[j][1]);      
//         //assert(mols[i].veloc[j] < MAX_VELOC);
      }   
      fixPosition(i);
     }
}

/**
 * Procedure: predictRotationsGear
 * Purpose: perform gear predictor step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void predictRotationsGear(void)
{
     int i, j;
     waterMol *m;
     double cr[] = {19.0, -10.0, 3.0};
     double cv[] = {27.0, -22.0, 7.0};
     double div = 24.0;

     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];

      // update quarternion 
      for (j = 0; j < DIM; j++)
      {
           m->qold[j] = m->q[j];
           m->qvold[j] = m->qv[j];

               // update positions, qvities, qaerations, and other time derivatives
           m->q[j] += props.deltaT*m->qv[j];
           m->q[j] += (props.deltaT2/div) * (cr[0]*m->qa[j] + cr[1]*m->qa1[j] + cr[2]*m->qa2[j]); 
           m->qv[j] = (m->q[j] - m->qold[j])/props.deltaT;
           m->qv[j] += (props.deltaT/div) * (cv[0]*m->qa[j] + cv[1]*m->qa1[j] + cv[2]*m->qa2[j]);
           m->qa2[j] = m->qa1[j];
           m->qa1[j] = m->qa[j];
      }
     }
}

/**
 * Procedure: correctRotationsGear
 * Purpose: perform gear corrector step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void correctRotationsGear(void)
{
     int i, j;
     waterMol *m;
     double cr[] = {3.0, 10.0, -1.0};
     double cv[] = {7.0, 6.0, -1.0};
     double div = 24.0;

     stats.AKE = 0;
     stats.AM[0] = stats.AM[1] = stats.AM[2] = 0;
     for (i = 0; i < numMols; i++)
     {
      m = &mols[i];

      // update velocities
      for (j = 0; j < DIM; j++)
      {
           m->q[j] = m->qold[j] + props.deltaT*m->qvold[j];
           m->q[j] += (props.deltaT2/div) * (cr[0]*m->qa[j] + cr[1]*m->qa1[j] + cr[2]*m->qa2[j]); 
           m->qv[j] = (m->q[j] - m->qold[j])/props.deltaT;
           m->qv[j] += (props.deltaT/div) * (cv[0]*m->qa[j] + cv[1]*m->qa1[j] + cv[2]*m->qa2[j]);

//         computeAngularVelocity(m);
           stats.AKE += m->inert[j]*m->omega[j]*m->omega[j];    // angular kinetic energy
           stats.AM[j] += m->omega[j];                          // angular momentum
      }   

      normalizeQuarternion(m->q);
     }
     stats.AKE *= 0.5;
}

/**
 * Procedure: correctRotationsVerlet
 * Purpose: perform verlet predictor step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void predictRotationsVerlet(void)
{
     int i, j;
     for (i = 0; i < numMols; i++)
     {
      // rotate sites
      setSites(&mols[i]);

      // update quarternion 
      for (j = 0; j < DIM; j++)
      {
           mols[i].q[j] += mols[i].qv[j]*props.deltaT + 0.5*mols[i].qa[j]*props.deltaT2;      
           mols[i].qv[j] += 0.5*mols[i].qa[j]*props.deltaT;
      }

      for (j = 0; j < NUM_SITES; j++)
           rotateToSpaceFrame(mols[i].sites[j].pos, mols[i].q);
     }
}

/**
 * Procedure: correctRotationsVerlet
 * Purpose: perform verlet corrector step for computing 
 *			Newton's differential equations of motion
 * Inputs: None
 * Outputs: None
 **/
void correctRotationsVerlet(void)
{
     int i, j;

     stats.AKE = 0;
     stats.AM[0] = stats.AM[1] = stats.AM[2] = 0;
     for (i = 0; i < numMols; i++)
     {
      // update velocities
      for (j = 0; j < DIM; j++)
      {
           mols[i].qv[j] += 0.5*mols[i].qa[j]*props.deltaT; // correct angular velocity

           stats.AKE += mols[i].inert[j]*mols[i].omega[j]*mols[i].omega[j];    // angular kinetic energy
           stats.AM[j] += mols[i].omega[j];                                       // angular momentum
      }   

      normalizeQuarternion(mols[i].q);
     }
     stats.AKE *= 0.5;
}

/**
 * Procedure: rotateMolecule
 * Purpose: rotate a molecule
 * Inputs: the molecule and the euler angle to rotate it by
 * Outputs: None
 **/
void rotateMolecule(waterMol *mol, double euler[DIM])
{
     int i;
     for (i = 0; i < NUM_SITES; i++)
      rotateEuler(euler, mol->sites[i].pos);
}

/**
 * Procedure: distance
 * Purpose: compute distance between two molecules
 * Inputs: the molecules and the distance
 * Outputs: the distance squared
 **/
double distance (waterMol *m1, waterMol *m2, double *dr)
{
     int k;
     double d2 = 0.0;
     for (k = 0; k < DIM; k++)
     {    
      dr[k] = m1->cm[k] - m2->cm[k];

      // periodic
      if (dr[k] > 0.5*bounds[k][1])
           dr[k] -= bounds[k][1];
      else if (dr[k] < -0.5*bounds[k][1])
           dr[k] += bounds[k][1];       

      d2 += dr[k] * dr[k];
     }
     return d2;
}

/**
 * Procedure: siteDistance
 * Purpose: compute distance between two interaction sites
 * Inputs: the molecules, the interaction sites, and the distance
 * Outputs: the distance squared
 **/
double siteDistance (waterMol *m1, waterMol *m2, site *s1, site *s2, double *dr)
{
     int k;
     double d2 = 0.0;
     for (k = 0; k < DIM; k++)
     {
      dr[k] = m1->cm[k] - m2->cm[k];

      // periodic
      if (dr[k] > 0.5*bounds[k][1])
           dr[k] -= bounds[k][1];
      else if (dr[k] < -0.5*bounds[k][1])
           dr[k] += bounds[k][1];       

      dr[k] += s1->pos[k] - s2->pos[k]; // add site distance

      d2 += dr[k] * dr[k];
     }
     return d2;
}

/**
 * Procedure: computeAlphas
 * Purpose: compute angular acceleration of a molecule (r x F)
 * Inputs: the molecule
 * Outputs: None
 **/
void computeAlphas(waterMol *m)
{
     double inert[DIM];
     int i, j, k;
     site *s;

     for (i = 0; i < DIM; i++)
      m->inert[i] = 0.0;

     for (i = 0; i < NUM_SITES; i++)
     {
      // use cross product to compute torque
      s = &m->sites[i];

      m->alpha[0] += s->pos[1] * s->force[2] - s->pos[2] * s->force[1];
      m->alpha[1] += s->pos[2] * s->force[0] - s->pos[0] * s->force[2];
      m->alpha[2] += s->pos[0] * s->force[1] - s->pos[1] * s->force[0];

      // compute moment of inertia
      for (i = 0; i < NUM_SITES; i++)
      {
           m->inert[0] += siteMass[i] * (s->pos[1] * s->pos[1] + s->pos[2] * s->pos[2]);
           m->inert[1] += siteMass[i] * (s->pos[0] * s->pos[0] + s->pos[2] * s->pos[2]);
           m->inert[2] += siteMass[i] * (s->pos[0] * s->pos[0] + s->pos[1] * s->pos[1]);
      }
     }
 
     for (i = 0; i < DIM; i++)
      m->alpha[i] /= m->inert[i];
}

/**
 * Procedure: computeSiteForces
 * Purpose: compute forces between two molecules
 * Inputs: the molecules, the potential energy, the kinetic energy
 * Outputs: None
 **/
void computeSiteForces(waterMol *m1, waterMol *m2, double *U, double *V)
{
     double f, df, u, dr[DIM];
     int i, j, k;
     site *s1, *s2;

     *U = *V = 0.0;

     // compute site forces
     for (i = 0; i < NUM_SITES; i++)
     {
      for (j = 0; j < NUM_SITES; j++)
      {
           s1 = &m1->sites[i];
           s2 = &m2->sites[j];

           computePairForce(m1, m2, s1, s2, dr, &f, &u);
           *U += u;

           for (k = 0; k < DIM; k++)
           {
            df = f * dr[k];
            
            // site forces for rotational motion
            s1->force[k] += df;
            s2->force[k] -= df;
      
            // total forces for translational motion
            m1->accel[k] += df;
            m2->accel[k] -= df;
            
//          *V += dr[k] * df;
            *V += df;
           }
      }
     }

//     printf("\n");
}

/**
 * Procedure: computePairForce
 * Purpose: compute the forces on a pair of interaction sites
 * Inputs: the molecules, the sites, the distance, the force, the potential 
 * Outputs: None
 **/
void computePairForce(waterMol *m1, waterMol *m2, site *s1, site *s2, double dr[DIM], double *f, double *u)
{
     int type;
     double r2, rm2, rm6;

     *f = *u = 0.0;

     type = s1->type | s2->type;

     r2 = siteDistance(m1, m2, s1, s2, dr);
     rm2 = 1.0/r2;

     switch (type)
     {
     case O_SITE:
      rm6 = rm2*rm2*rm2;
      *f = 48.0 * rm6 * ( rm6 - 0.5 ) * rm2; 
      *u = 4.0 * rm6 * ( rm6 - 1.0 );
//    *f = 48.0 * rm6 * ( LJ_A * rm6 - 0.5 * LJ_C ) * rm2; 
//    *u = 4.0 * rm6 * ( LJ_A * rm6 - LJ_C );
//    printf("OO ");
      break;

#if COULOMB_POTENTIALS

     case M_SITE:
      *u = 4.0 * B_COULOMB * sqrt (rm2);
//    *u = C_QM * C_QM * sqrt (rm2);
      *f += *u * rm2;
//    printf("MM ");
      break;

     case H_SITE:
      *u = B_COULOMB * sqrt (rm2);
//    *u = C_QH * C_QH * sqrt (rm2);
      *f += *u * rm2;
//    printf("HH ");
      break;

     case MH_SITE:
      *u = -2.0 * B_COULOMB * sqrt (rm2);
//    *u = C_QM * C_QH * sqrt (rm2);
      *f += *u * rm2;
//    printf("MH ");
      break;

#endif

     default:
//    printf("XX ");
      break;
     }
}
