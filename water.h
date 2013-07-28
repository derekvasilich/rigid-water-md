/**
 * Module: md.c
 * Author: Derek Williams
 * Created: 01/30/04
 * Purpose: 3-D simulation of liquid water using the TIP4P model
 */

#ifndef WATER_H
#define WATER_H

#define ARGON                  0
#define ALL_PAIRS              1
#define ADJUST_TEMP_STEP     100
#define EQUILIBERATION_TIME  4000
#define NUM_MOLECULES        32//108//864//108//864//108//108//864//2048//108//864//32//108//2048
#define DENSITY                0.98//0.844//0.98
#define TEMPERATURE            3.8//0.71//3.8

#if ARGON
    #define USE_GEAR            0
    #define ROTATIONAL_MOTION   0
    #define COULOMB_POTENTIALS  0
    #define CUTOFF_DISTANCE     2.8         // maximum distance to compute potential
    #define DT                  0.005       // time step
#else
    #define USE_GEAR            0
    #define ROTATIONAL_MOTION   1
    #define COULOMB_POTENTIALS  1
    #define CUTOFF_DISTANCE     100000//2.8         // maximum distance to compute potential
    #define DT                  0.0005       // time step
#endif

#define RDF_BINS 100

#define DATA_FILE_NAME      "water.bin"

#define BALL_AND_STICK      0
#define FIXED_SITES         0
#define MAX_VELOC           100.0

#define SIZE_X              600          // width of window in pixels
#define SIZE_Y              600          // height of window in pixels

#define OXYGEN_MASS         8.0/9.0      // mass of oxygen atom in reduced units
#define HYDROGEN_MASS       1.0/18.0     // mass of hydrogen atom in reduced 

#define OXYGEN_RADIUS       SPHERE_RADIUS
#define HYDROGEN_RADIUS     cbrt(HYDROGEN_MASS)*OXYGEN_RADIUS
#define B_COULOMB           183.5        // should be 183.5

// DO NOT CHANGE
#define DIM                 3            // dimension of system
////////////////

#define NUM_SITES           4

#if BALL_AND_STICK
#define SPHERE_RADIUS       0.15
#else
#define SPHERE_RADIUS       0.3
#endif

#define COS(x) cos(x*M_PI/180.0)
#define SIN(x) sin(x*M_PI/180.0)

#define SQR(X) X*X

#ifndef cbrt
#define cbrt(x)   pow(x, 1.0/3.0)
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// 0.5 * (quarternion rotation matrix)
#define A11(q) SQR(q[0]) + SQR(q[3]) - 0.5 
#define A12(q) q[0]*q[1] + q[2]*q[3]       
#define A13(q) q[0]*q[2] - q[1]*q[3]
#define A21(q) q[0]*q[1] - q[2]*q[3]                   
#define A22(q) SQR(q[1]) + SQR(q[3]) - 0.5
#define A23(q) q[1]*q[2] + q[0]*q[3]
#define A31(q) q[0]*q[2] + q[1]*q[3]                   
#define A32(q) q[1]*q[2] - q[0]*q[3]       
#define A33(q) SQR(q[2]) + SQR(q[3]) - 0.5

// W matrix for computing quaternion accelerations
#define W11(Q) Q[3]
#define W12(Q) Q[2]  
#define W13(Q) -Q[1]  
#define W14(Q) -Q[0]
#define W21(Q) -Q[2]  
#define W22(Q) Q[3]  
#define W23(Q) Q[0]  
#define W24(Q) Q[1]
#define W31(Q) Q[1]  
#define W32(Q) -Q[0]  
#define W33(Q) Q[3]  
#define W34(Q) -Q[2]
#define W41(Q) Q[0]  
#define W42(Q) Q[1]  
#define W43(Q) Q[2]  
#define W44(Q) Q[3]

struct site_t
{
     int type;                            // type of site
     double force[DIM];                   // force on site
     double pos[DIM];                     // position of site relative to center of mass of molecule
};

typedef struct site_t site;

struct waterMol_t
{
     int n;
     site sites[NUM_SITES];                // the sites
     int cell[3];                          // current cell assignment
     nodeptr node;                         // ptr to cell assignment node
     double cm[DIM];                       // position of the center of mass of the molecule
     double cmold[DIM];                    // position of the center of mass of the molecule
     double inert[DIM];                    // moment of inertia
     double euler[DIM];                    // euler angles specifying the orientation of the molecules
     double alpha[DIM];                    // angular acceleration
     double omega[DIM];                    // angular velocity
     double accel[DIM];                    // linear acceleration
     double accel1[DIM];                   // 1st derivative of linear acceleration
     double accel2[DIM];                   // 2nd derivative of linear acceleration
     double veloc[DIM];                    // linear velocity
     double velocold[DIM];                 // linear velocity
     double tau[DIM];                      // torque
     double q[4];                          // quarternion orientation
     double qv[4];                         // quarternion velocity
     double qold[4];                       // quarternion orientation
     double qvold[4];                      // quarternion velocity
     double qa[4];                         // quarternion acceleration
     double qa1[4];                        // 1st derivative of quarternion acceleration
     double qa2[4];                        // 2nd derivative of quarternion acceleration
};

typedef struct waterMol_t waterMol;

double siteMass[NUM_SITES];

int numMols;                        // number of molecules
waterMol *mols;                     // the molecules
double *distances;
int pairs;

nodeptr ***cells;
int numCells[3];
double cellSize;

#define M_SITE         0x0001       // identify M-sites
#define O_SITE         0x0010       // identify oxygen sites
#define H_SITE         0x0100       // identify hydrogen sites
#define MH_SITE        0x0101       // identify hydrogen M-site pair

int width;
int height;

double bounds[3][2];

struct stat_t
{
     double AKE;                     // instantanious angular kinetic energy
     double KE;                      // instantanious kinetic energy
     double UE;                      // instantanious potential energy
     double LM[DIM];                 // linear momentum of the system
     double AM[DIM];                 // angular momentum of the system
     double VIRIAL;                  // virial
     double accVIRIAL;               // virial accumulator
     double accKE;                   // kinetic energy accumulator
     double accAKE;                  // rotation kinetic energy accumulator
     double T;                       // temperature
     double P;                       // presure
     double V;                       // volume
     double meanKE;                  // mean kinetic energy
     double meanAKE;                 // mean rotational kinetic energy
     int accRDFOO[RDF_BINS];                       // accumulator for RDF
     int accRDFHH[RDF_BINS];                       // accumulator for RDF
     int accRDFOH[RDF_BINS];                       // accumulator for RDF
     double nRDF;
     double rdfMax;
     double dr;
     double accN;
     double time;
     int step;
};

typedef struct stat_t stat;

stat stats;

struct prop_t
{
     double rCut;
     double rrCut;
     int useGL;                          // render?
     int displayStep;                    // number of steps to display
     int tempAdjustStep;
     int maxTime;                        // maximum number of time steps
     double deltaT;
     double deltaT2;
     double vMag;
     double vqMag;
     int hydrogenSphere;
     int oxygenSphere;
     double cutoff2;
     double sigma;
     double epsilon;
     double density;
     int equiliberationSteps;
};

typedef struct prop_t prop;

prop props;

#define ENERGY_STATS                1
#define VELOCITY_X_STATS            2
#define VELOCITY_Y_STATS            3
#define VELOCITY_Z_STATS            4
#define ANGULAR_VELOCITY_STATS      5
#define FPS_STATS                   6
#define MOMENTUM_STATS              7
#define MOLECULE_STATS              8
#define PRESSURE_TEMPERATURE_STATS  9
#define PAIR_DISTRIBUTION_STATS    10
#define QUATERNION_STATS           11
#define RDF_STATS                  12

#define ENERGY_FILE			"energy.dat"
#define VELOCITY_FILE		"velocity.dat"
#define ANGULAR_FILE		"angularvelocity.dat"
#define MOMENTUM_FILE		"momentum.dat"
#define PRESSURETEMP_FILE	"pressuretemp.dat"
#define RDF_FILE			"rdf.dat"

#define USE_GL_DEFAULT				1
#define DISPLAY_STEP_DEFAULT		10
#define MAX_TIME_DEFAULT			10000

struct statfile_t
{
     FILE *energy;
     FILE *velocity;  
     FILE *angular;
     FILE *momentum;
     FILE *pressuretemp;
     FILE *rdf;
};

typedef struct statfile_t statfile;

statfile statfiles;

int statsType;

char* helpStr = "water dispStep totalStep [statsType] [useGL]\n\nDescription:\n\tSimulate liquid water using a 4 site rigid molecular model.\n\nOptions:\n\tdispStep - the number of time steps until redisplay\n\ttotalStep - the total number of time steps for the simulation\n\tstatsType (optional) - the type of stats to display.\n\t\t 1 - Energy\n\t\t 2 - Center of Mass Velocities\n\t\t 4 - Angular Velocities\n\t\t 8 - Frames Per Second\n\t\t16 - Momentum\n\tuseGL (optional) - 0 or 1 indicating weather or not to render using openGL\n";

time_t start;

int main(int argc, char *argv[]);
int processArguments(int argc, char *argv[]);

int save(char *fileName);
int load(char *fileName);

void assignCells(void);
void reassignCells(void);
void initializeCells(void);

void resetAccumulators(void);

static void draw(void);
static void reshape(int w, int h);
void initializeGL(void);

static void idle(void);
void predictVerlet(void);
void correctVerlet(void);
void predictGear(void);
void correctGear(void);
void predictRotationsVerlet(void);
void correctRotationsVerlet(void);
void predictRotationsGear(void);
void correctRotationsGear(void);
int initialize(void);
void buildLattice(void);
void intializeVelocities(void);
void adjustTemperature(void);
void computeInertia(waterMol *mol);
void initializeCells(void);
void drawMolecule(waterMol* mol);
void evaluateForces(void);
void computeRotationalAccelerations(void);
void computeSiteForces(waterMol *m1, waterMol *m2, double *U, double *V);
void computePairForce(waterMol *m1, waterMol *m2, site *s1, site *s2, double dr[DIM], double *f, double *u);
void computeAlphas(waterMol *m);
void computeTorque(waterMol *m);
void rotateToBodyFrame(double v[DIM], double q[4]);
void rotateToSpaceFrame(double v[DIM], double q[4]);
void computeAngularVelocity(waterMol *m);
void randomQuarternion(double q[4]);
void normalizeQuarternion(double q[4]);
double distance(waterMol *m1, waterMol *m2, double *dr);
double siteDistance (waterMol *m1, waterMol *m2, site *s1, site *s2, double *dr);
void drawBoundary(void);
void fixPosition(int i);
void displayStats(FILE* fp, int type);
void initializeMolecule(waterMol *mol);
void setSites(waterMol *mol);
void randomUnitVector(double vec[DIM]);
void rotateEuler(double euler[DIM], double v[DIM]);
void rotateMolecule(waterMol *mol, double euler[DIM]);
void updateStats(void);
void initializePropsStats(void);
void openStatsFiles(void);
void closeStatsFiles(void);
void doStats(void);
void computeKE(void);
void initializeAngularVelocities(void);
void predict(void);
void correct(void);
void predictRotations(void);
void correctRotations(void);

#endif // WATER_H
