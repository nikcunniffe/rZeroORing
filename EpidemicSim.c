#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <direct.h>
#include <process.h>

/* MT19937 random number generation */
#include "mt19937ar.h"

/* allow Visual Studio to compile ANSI C containing strcpy */
#pragma warning(disable : 4996)

#define	_MAX_STR_LEN		1024
#define	_STR_SEP			",\t "
#define _BLOCK_SIZE			256
#define	_PI					3.1415926535897932384626433
#define	_NOT_SET			-1
#define _NUMERIC_UNDERFLOW	1e-5
#define N_DUMP_STEPS		100			/* used if dumping out time courses rather than generations */
#define	_ONE_LINE_GEN_OUT	1			/* whether or not to put all information for a generation on a single line */

#ifdef _WIN32
#define 	C_DIR_DELIMITER '\\'
#else
#define 	C_DIR_DELIMITER '/'
#endif

enum
{
	SUSCEPTIBLE = 0,
	INFECTED = 1,
	REMOVED = 2
} hostStatus;

enum
{
	TYPE_I = 1,
	TYPE_II = 2
} hostType;

enum
{
	KERNEL_I = 1,		/* exponential power */
	KERNEL_II = 2		/* flat (for testing) */
} kernelType;

enum
{
	DUMP_GENS = 1,
	DUMP_TIMES = 2
} dumpType;

enum
{
	MODEL_SIS = 1,
	MODEL_SIR = 2
} modelType;

typedef struct {
	double	dThetaOne;		/* Infectivity */
	double	dThetaTwo;
	double	dRhoOne;		/* Susceptibility */
	double	dRhoTwo;
	double	dMuOne;			/* Death rate */
	double	dMuTwo;
	int		nInitOne;		/* Numbers of initial infections */
	int		nInitTwo;
	int		bCacheKernel;	/* Whether or not to store the kernel in memory */
	int		eKernelType;	/* What sort of kernel */
	double	dA;				/* Exponential-power kernel */
	double	dC;
	int		nNumIts;		/* Number of iterations to run */
	int		nMaxGen;		/* Make infection rate = 0 after this many generations */
	int		eModelType;		/* Whether SIS or SIR model */
	char	sXYFile[_MAX_STR_LEN];
	char	sOutFile[_MAX_STR_LEN];
	char	sParamDumpFile[_MAX_STR_LEN];
	double	dMaxTime;
	int		eDumpType;
	int		bDumpHostStatus;
} t_Params;

typedef struct {
	double	dX;
	double	dY;
	int		eType;
} t_SingleHost;

typedef struct {
	t_SingleHost	*aHosts;
	int				nHosts;
	int				nAlloc;
} t_Hosts;

typedef struct {
	int		eStatus;
	double	dRate;
	int		nGen;
	int		nEntryPtr;	/*  only for infected hosts, store where the
							relevant entry is in the epidemicEntryList */
} t_HostStatus;

typedef struct {
	int		nGen;
	double	dInfectTime;
	double	dRemovalTime;
	int		nHostID;
	int		eType;
} t_EpidemicEntry;

typedef struct {
	t_EpidemicEntry *aEntries;
	int				nEntries;
	int				nAlloc;
} t_Epidemic;

typedef struct {
	double	*aKernel;	/* stored as a flattened array */
} t_Kernel;

/*
	seed random number generator
*/
void	seedRandom(unsigned long ulnSeed)
{
	unsigned long		myPID;

	if (ulnSeed <= 0)
	{
		ulnSeed = (unsigned long)time(NULL);
#ifndef _WIN32	/* make sure different processes started at same time have different seeds */
		myPID = (unsigned long)getpid();
#else
		myPID = (unsigned long)_getpid();
#endif
		ulnSeed += myPID;
	}
	init_genrand(ulnSeed);
}

/* work out configuration file name and check whether it exists */
int	getCfgFileName(char *szProgName, char *szCfgFile)
{
	char	*pPtr;
	FILE 	*fp;
	char	szDir[_MAX_STR_LEN];

	szCfgFile[0] = '\0';
	{
		if ((pPtr = strrchr(szProgName, C_DIR_DELIMITER)) != NULL)
		{
			strcpy(szCfgFile, pPtr + 1);
		}
		else
		{
			strcpy(szCfgFile, szProgName);
		}
		if ((pPtr = strstr(szCfgFile, ".exe")) != NULL)
		{
			*pPtr = '\0';
		}
		strcat(szCfgFile, ".cfg");
	}
	/* check file exists */
	fp = fopen(szCfgFile, "rb");
	if (fp)
	{
		fclose(fp);
		return 1;
	}
	getcwd(szDir, _MAX_STR_LEN);
	fprintf(stderr, "Did not find config file %s in %s\n", szCfgFile, szDir);
	return 0;
}

/*
(rather slow)

routines to find values from the command line options,
or, failing that, from the cfg file
*/
int findKey(int argc, char **argv, char*szCfgFile, char *szKey, char *szValue)
{
	char *pVal;
	int	 bRet, i;
	FILE *fp;
	char *pThisPair;
	char *szArgvCopy;

	bRet = 0;
	i = 0;
	/* try to find the relevant key on the command line */
	while (bRet == 0 && i<argc)
	{
		szArgvCopy = strdup(argv[i]);
		if (szArgvCopy)
		{
			pThisPair = strtok(szArgvCopy, " \t");
			while (pThisPair)
			{
				if (strncmp(pThisPair, szKey, strlen(szKey)) == 0)
				{
					pVal = strchr(pThisPair, '=');
					if (pVal)
					{
						if (pThisPair[strlen(szKey)] == '=') /* check full string matches the key */
						{
							strcpy(szValue, pVal + 1);
							fprintf(stdout, "extracted %s->%s from command line\n", szKey, szValue);
							bRet = 1;
						}
					}
				}
				pThisPair = strtok(NULL, " \t");
			}
			free(szArgvCopy);
		}
		i++;
	}
	/* otherwise, look in the cfg file */
	if (bRet == 0)
	{
		fp = fopen(szCfgFile, "rb");
		if (fp)
		{
			char szLine[_MAX_STR_LEN];

			while (!bRet && fgets(szLine, _MAX_STR_LEN, fp))
			{
				char *pPtr;
				if ((pPtr = strchr(szLine, '=')) != NULL)
				{
					*pPtr = '\0';
					if (strcmp(szKey, szLine) == 0)
					{
						strcpy(szValue, pPtr + 1);
						/* strip off newline (if any) */
						if ((pPtr = strpbrk(szValue, "\r\n")) != NULL)
							*pPtr = '\0';
						bRet = 1;
					}
				}
			}
			fclose(fp);
		}
	}
	return bRet;
}

int readStringFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, char *szValue)
{
	return(findKey(argc, argv, szCfgFile, szKey, szValue));
}

int readDoubleFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, double *pdValue)
{
	char szValue[_MAX_STR_LEN];

	if (findKey(argc, argv, szCfgFile, szKey, szValue))
	{
		*pdValue = atof(szValue);
		return 1;
	}
	return 0;
}

int readIntFromCfg(int argc, char **argv, char *szCfgFile, char *szKey, int *pnValue)
{
	char szValue[_MAX_STR_LEN];

	if (findKey(argc, argv, szCfgFile, szKey, szValue))
	{
		*pnValue = atoi(szValue);
		return 1;
	}
	return 0;
}

/*
	return uniform number on [0,1)
*/
double	uniformRandom()
{
	return genrand_real3();
}

/*
	position in the flattened array for a pair of hosts
*/
int	posFromHostIDs(int idOne, int idTwo, int numHosts)
{
	return idOne + numHosts*idTwo;
}

/*
	Logarithmic gamma function using the algorithm from Numerical Recipes
*/ 
double logGamma(double x) 
{

	int				i;
	double			a, b, t, s;
	static double	c[6] = { 76.18009172947146,
							-86.50532032941677,
							24.01409824083091,
							-1.231739572450155,
							0.1208650973866179e-2,
							-0.5395239384953e-5};
	
	b = a = x;
	t = a + 5.5;
	t -= (a + 0.5) * log(t);
	s = 1.000000000190015;
	for (i = 0; i <= 5; i++)
	{
		s += c[i] / ++b;
	}

	return -t + log(2.5066282746310005 * s / a);
}

double simpleGammaFunction(double x)
{
	double retVal;
	
	retVal = exp(logGamma(x));
	return retVal;
}

/*
	dispersal kernel
*/
double dispKernel(double r, double alpha, double c, int eKernelType)
{
	double dK;
	double dExp;

	dK = 1.0;
	if (eKernelType == KERNEL_I)
	{
		dExp = pow(r / alpha, c);
		dK = c * exp(-dExp) / (2.0 * _PI * alpha * alpha * simpleGammaFunction(2.0 / c));
	}
	return dK;
}

/*
	read in parameters (will alter to read in from cfg file)
*/
int setParams(t_Params *pParams)
{
	pParams->dThetaOne = 0.1;			/* infectivity of a type 1 host */
	pParams->dThetaTwo = 0.05;			/* infectivity of a type 2 host */
	pParams->dRhoOne = 0.05;			/* susceptibility of a type 1 host  */
	pParams->dRhoTwo = 0.05;			/* susceptibility of a type 2 host  */
	pParams->dMuOne = 2.5;				/* death rate of a type 1 host  */
	pParams->dMuTwo = 1.0;				/* death rate of a type 2 host  */
	pParams->nInitOne = 1;				/* initial number of type 1 infections */
	pParams->nInitTwo = 1;				/* initial number of type 2 infections */
	pParams->bCacheKernel = 1;			/* kernel parameters */
	pParams->eKernelType = KERNEL_I;	/* exponential power kernel */
	pParams->dA = 0.5;					/* scale */
	pParams->dC = 1.0;					/* power */
	pParams->nNumIts = 100;
	/* in this generation the infection rate is artificially set to zero */
	pParams->nMaxGen = 3;//5;
	strcpy(pParams->sXYFile, "finalExampleLS_xy.csv");
	strcpy(pParams->sOutFile, "finalExampleLS_Epidemics_Eg5_1.csv");
	pParams->dMaxTime = 20.0;
	pParams->eDumpType = DUMP_GENS;
	return 1;
}

int dumpParametersToCSV(t_Params *pParams)
{
	FILE *fOut;

	fOut = fopen(pParams->sParamDumpFile, "wb");
	if (!fOut)
	{
		fprintf(stderr, "dumpParametersToCSV(): could not open file\n");
		return 0;
	}
	fprintf(fOut, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
		"thetaOne",
		"thetaTwo",
		"rhoOne",
		"rhoTwo",
		"muOne",
		"muTwo",
		"initOne",
		"initTwo",
		"kernelType",
		"dispA",
		"dispC",
		"numIts",
		"maxGen",
		"xyFile",
		"modelType");
	fprintf(fOut, "%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%d,%d,%d,%.7f,%.7f,%d,%d,%s,%d\n",
		pParams->dThetaOne,
		pParams->dThetaTwo,
		pParams->dRhoOne,
		pParams->dRhoTwo,
		pParams->dMuOne,
		pParams->dMuTwo,
		pParams->nInitOne,
		pParams->nInitTwo,
		pParams->eKernelType,
		pParams->dA,
		pParams->dC,
		pParams->nNumIts,
		pParams->nMaxGen,
		pParams->sXYFile,
		pParams->eModelType);
	fclose(fOut);
	return 1;
}

int readParams(t_Params *pParams, int argc, char **argv)
{
	char szCfgFile[_MAX_STR_LEN];

	fprintf(stdout, "readParams()\n");
	memset(pParams, 0, sizeof(t_Params));
	pParams->bCacheKernel = 1;			/* kernel parameters */
	pParams->eDumpType = DUMP_GENS;		/* dump out generations only */
	pParams->dMaxTime = -1;
	if (!getCfgFileName(argv[0], szCfgFile))
	{
		fprintf(stderr, "readParams(): Couldn't find cfg file for program name '%s'\n", argv[0]);
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "thetaOne", &pParams->dThetaOne))
	{
		fprintf(stderr, "readParams(): Couldn't read thetaOne\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "thetaTwo", &pParams->dThetaTwo))
	{
		fprintf(stderr, "readParams(): Couldn't read thetaTwo\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "rhoOne", &pParams->dRhoOne))
	{
		fprintf(stderr, "readParams(): Couldn't read rhoOne\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "rhoTwo", &pParams->dRhoTwo))
	{
		fprintf(stderr, "readParams(): Couldn't read rhoTwo\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "muOne", &pParams->dMuOne))
	{
		fprintf(stderr, "readParams(): Couldn't read muOne\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "muTwo", &pParams->dMuTwo))
	{
		fprintf(stderr, "readParams(): Couldn't read muTwo\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "initOne", &pParams->nInitOne))
	{
		fprintf(stderr, "readParams(): Couldn't read initOne\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "initTwo", &pParams->nInitTwo))
	{
		fprintf(stderr, "readParams(): Couldn't read initTwo\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "kernelType", &pParams->eKernelType))
	{
		fprintf(stderr, "readParams(): Couldn't read kernelType\n");
		return 0;
	}
	if (!(pParams->eKernelType == KERNEL_I || pParams->eKernelType == KERNEL_II))
	{
		fprintf(stderr, "readParams(): Invalid kernelType (must be %d or %d)\n", KERNEL_I, KERNEL_II);
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "dispA", &pParams->dA))
	{
		fprintf(stderr, "readParams(): Couldn't read dispA\n");
		return 0;
	}
	if (!readDoubleFromCfg(argc, argv, szCfgFile, "dispC", &pParams->dC))
	{
		fprintf(stderr, "readParams(): Couldn't read dispC\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "numIts", &pParams->nNumIts))
	{
		fprintf(stderr, "readParams(): Couldn't read numIts\n");
		return 0;
	}
	if (!readIntFromCfg(argc, argv, szCfgFile, "maxGen", &pParams->nMaxGen))
	{
		fprintf(stdout, "readParams(): Couldn't read maxGen\n");
		return 0;
	}
	if (!readStringFromCfg(argc, argv, szCfgFile, "xyFile", pParams->sXYFile))
	{
		fprintf(stdout, "readParams(): Couldn't read xyFile\n");
		return 0;
	}
	if (!readStringFromCfg(argc, argv, szCfgFile, "outFile", pParams->sOutFile))
	{
		fprintf(stderr, "readParams(): Couldn't read outFile\n");
		return 0;
	}
	/* model type: SIS (default) or SIR: this field is not required in config file */
	pParams->eModelType = MODEL_SIS;
	readIntFromCfg(argc, argv, szCfgFile, "modelType", &pParams->eModelType);
	if (!(pParams->eModelType == MODEL_SIS || pParams->eModelType == MODEL_SIR))
	{
		fprintf(stderr, "readParams(): Invalid modelType (must be %d or %d)\n", MODEL_SIS, MODEL_SIR);
		return 0;
	}
	/* whether or not to dump information on host status...note is not required */
	pParams->bDumpHostStatus = 0;
	readIntFromCfg(argc, argv, szCfgFile, "dumpHostStatus", &pParams->bDumpHostStatus);
	/* create filename for dump of all parameters and actually do the dump */
	{
		char *sTmp, *p;

		sTmp = strdup(pParams->sOutFile);
		if (!sTmp)
		{
			fprintf(stderr, "readParams(): Out of memory\n");
			return 0;
		}
		/* strip any extension */
		if (p = strrchr(sTmp, '.'))
		{
			*p = '\0';
		}
		sprintf(pParams->sParamDumpFile, "%s_param.csv", sTmp);
		free(sTmp);
	}
	return dumpParametersToCSV(pParams);
}

/*
	calculate and store the dispersal kernel
*/
int calcKernel(t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel)
{
	int		p,i,j,retVal;
	double	d,k;

	retVal = 1;
	if (pParams->bCacheKernel)
	{
		retVal = 0;
		pKernel->aKernel = malloc(sizeof(double)*pHosts->nHosts*pHosts->nHosts);
		if (pKernel->aKernel)
		{
			/*
				set kernel between pairs of hosts
				note copies top triangle to the bottom rather than recalculating
			*/
			for (i = 0; i < pHosts->nHosts; i++)
			{
				for (j = 0; j < i; j++)
				{
					d = pow(pHosts->aHosts[i].dX - pHosts->aHosts[j].dX, 2);
					d += pow(pHosts->aHosts[i].dY - pHosts->aHosts[j].dY, 2);
					d = sqrt(d);
					k = dispKernel(d, pParams->dA, pParams->dC, pParams->eKernelType);
					p = posFromHostIDs(i, j, pHosts->nHosts);
					pKernel->aKernel[p] = k;
					p = posFromHostIDs(j, i, pHosts->nHosts);
					pKernel->aKernel[p] = k;
					//fprintf(stdout, "%d %d %d %.5f %.5f\n", i, j, p, d, k);
				}
			}
			/* set kernel from a single host onto itself to be zero */
			for (i = 0; i < pHosts->nHosts; i++)
			{
				p = posFromHostIDs(i, i, pHosts->nHosts);
				pKernel->aKernel[p] = 0.0;
			}
			fprintf(stdout, "Set up kernel\n");
			retVal = 1;
		}
	}
	else
	{
		pKernel->aKernel = NULL;
	}
	return retVal;
}

/*
	read in host locations and types
*/
int loadHosts(t_Params *pParams, t_Hosts *pHosts)
{
	int		retVal,tokNum;
	FILE	*f;
	char	sBuff[_MAX_STR_LEN];
	char	*p;

	retVal = 0;
	pHosts->nHosts = 0;
	f = fopen(pParams->sXYFile, "rb");
	if (f)
	{
		/* discard the first line, which is just a header */
		fgets(sBuff, _MAX_STR_LEN, f);
		retVal = 1;
		/* loop around reading and parsing the rest of the file */
		while (retVal && fgets(sBuff, _MAX_STR_LEN, f))
		{
			if (pHosts->nAlloc == pHosts->nHosts)
			{
				pHosts->nAlloc += _BLOCK_SIZE;
				pHosts->aHosts = realloc(pHosts->aHosts, sizeof(t_SingleHost)* pHosts->nAlloc);
				if (!pHosts->aHosts)
				{
					retVal = 0;
				}
			}
			if (retVal)
			{
				tokNum = 0;
				p = strtok(sBuff, _STR_SEP);
				while (p)
				{
					switch (tokNum)
					{
					case 0:
						pHosts->aHosts[pHosts->nHosts].dX = atof(p);
						break;
					case 1:
						pHosts->aHosts[pHosts->nHosts].dY = atof(p);
						break;
					case 2:
						pHosts->aHosts[pHosts->nHosts].eType = atoi(p);
						if (!(pHosts->aHosts[pHosts->nHosts].eType == TYPE_I
								|| pHosts->aHosts[pHosts->nHosts].eType == TYPE_II))
						{
							/* invalid host type */
							retVal = 0;
						}
						pHosts->nHosts++;
						break;
					default:
						retVal = 0;
						break;
					}
					tokNum++;
					p = strtok(NULL, _STR_SEP);
				}
				if(tokNum != 3) /* after the while loop, tokNum should be one larger than the number of tokens per line */
				{
					retVal = 0;
				}
			}
		}
		fclose(f);
	}
	fprintf(stdout, "Read in %d hosts\n", pHosts->nHosts);
	return (pHosts->nHosts > 0);
}

double getKernel(int hostOne, int hostTwo, t_Kernel *pKernel, t_Hosts *pHosts, t_Params *pParams)
{
	int		p;
	double	dKernel;

	if (pParams->bCacheKernel)
	{
		p = posFromHostIDs(hostOne, hostTwo, pHosts->nHosts);
		dKernel = pKernel->aKernel[p];
	}
	else
	{
		fprintf(stderr, "Not yet implemented!\n");
		dKernel = _NOT_SET;
	}
	return dKernel;
}

int recoverHost(int thisHost, double thisTime, t_Epidemic *pEpidemic, t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel, t_HostStatus *aHostStatus, double *pTotalRate)
{
	int		i, retVal;
	double	thisTheta, thisRho, thisExtra;

	retVal = 1;
	if (pHosts->aHosts[thisHost].eType == TYPE_I)
	{
		*pTotalRate -= pParams->dMuOne;
		thisTheta = pParams->dThetaOne;
	}
	else
	{
		*pTotalRate -= pParams->dMuTwo;
		thisTheta = pParams->dThetaTwo;
	}
	/* update susceptible hosts to no longer feel the force of infection from this one */
	/* note only need to do this when host isn't so old that not infecting anyway */
	if (aHostStatus[thisHost].nGen < pParams->nMaxGen)
	{
		for (i = 0; i < pHosts->nHosts; i++)
		{
			if (aHostStatus[i].eStatus == SUSCEPTIBLE)
			{
				thisRho = pParams->dRhoOne;
				if (pHosts->aHosts[i].eType == TYPE_II)
				{
					thisRho = pParams->dRhoTwo;
				}
				thisExtra = thisTheta * thisRho * getKernel(i, thisHost, pKernel, pHosts, pParams);
				aHostStatus[i].dRate -= thisExtra;
				if (aHostStatus[i].dRate < 0.0)
				{
					aHostStatus[i].dRate = 0.0;
				}
				*pTotalRate -= thisExtra;
			}
		}
	}
	pEpidemic->aEntries[aHostStatus[thisHost].nEntryPtr].dRemovalTime = thisTime;

	if (pParams->eModelType == MODEL_SIS)
	{
		aHostStatus[thisHost].eStatus = SUSCEPTIBLE;
		aHostStatus[thisHost].dRate = 0.0;

		/* loop around and add force back onto this one from all infected hosts */
		thisRho = pParams->dRhoOne;
		if (pHosts->aHosts[thisHost].eType == TYPE_II)
		{
			thisRho = pParams->dRhoTwo;
		}
		for (i = 0; i < pHosts->nHosts; i++)
		{
			if (aHostStatus[i].eStatus == INFECTED && aHostStatus[i].nGen < pParams->nMaxGen)
			{
				thisTheta = pParams->dThetaOne;
				if (pHosts->aHosts[i].eType == TYPE_II)
				{
					thisTheta = pParams->dThetaTwo;
				}
				thisExtra = thisTheta * thisRho * getKernel(i, thisHost, pKernel, pHosts, pParams);
				aHostStatus[thisHost].dRate += thisExtra;
				*pTotalRate += thisExtra;
			}
		}
	}
	else
	{
		aHostStatus[thisHost].eStatus = REMOVED;
		aHostStatus[thisHost].dRate = 0.0;
	}
	if (*pTotalRate < 0.0)
	{
		*pTotalRate = 0.0;
	}
	return retVal;
}

int infectHost(int thisHost, double thisTime, int infectedBy, t_Epidemic *pEpidemic, t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel, t_HostStatus *aHostStatus, double *pTotalRate)
{
	int			retVal,i,thisGen;
	double		thisTheta,thisRho,thisExtra;

	retVal = 1;

	/* update this host's status */
	if (infectedBy >= 0)
	{
		thisGen = aHostStatus[infectedBy].nGen + 1;
	}
	else
	{
		thisGen = 0;
	}
	aHostStatus[thisHost].nGen = thisGen;
	*pTotalRate -= aHostStatus[thisHost].dRate;
	aHostStatus[thisHost].nGen = thisGen;
	if (pHosts->aHosts[thisHost].eType == TYPE_I)
	{
		aHostStatus[thisHost].dRate = pParams->dMuOne;
		thisTheta = pParams->dThetaOne;
	}
	else
	{
		aHostStatus[thisHost].dRate = pParams->dMuTwo;
		thisTheta = pParams->dThetaTwo;
	}
	if (thisGen >= pParams->nMaxGen)
	{
		thisTheta = 0.0;	/* artificially stop infections once too many generations have passed */
	}
	*pTotalRate += aHostStatus[thisHost].dRate;
	aHostStatus[thisHost].eStatus = INFECTED;
	/* update susceptible hosts to feel the new force of infection from this one */
	for (i = 0; i < pHosts->nHosts; i++)
	{
		if (aHostStatus[i].eStatus == SUSCEPTIBLE)
		{
			thisRho = pParams->dRhoOne;
			if (pHosts->aHosts[i].eType == TYPE_II)
			{
				thisRho = pParams->dRhoTwo;
			}
			thisExtra = thisTheta * thisRho * getKernel(i, thisHost, pKernel, pHosts, pParams);
			aHostStatus[i].dRate += thisExtra;
			*pTotalRate += thisExtra;
		}
	}
	/* update the epidemic information */
	if (pEpidemic->nAlloc == pEpidemic->nEntries)
	{
		pEpidemic->nAlloc += _BLOCK_SIZE;
		pEpidemic->aEntries = realloc(pEpidemic->aEntries, sizeof(t_EpidemicEntry)* pEpidemic->nAlloc);
		if (!pEpidemic->aEntries)
		{
			retVal = 0;
		}
	}
	if (retVal)
	{
		pEpidemic->aEntries[pEpidemic->nEntries].nGen = thisGen;
		pEpidemic->aEntries[pEpidemic->nEntries].dInfectTime = thisTime;
		pEpidemic->aEntries[pEpidemic->nEntries].eType = pHosts->aHosts[thisHost].eType;
		pEpidemic->aEntries[pEpidemic->nEntries].nHostID = thisHost;
		pEpidemic->aEntries[pEpidemic->nEntries].dRemovalTime = _NOT_SET;
		aHostStatus[thisHost].nEntryPtr = pEpidemic->nEntries;
		pEpidemic->nEntries++;
	}
	return retVal;
}

int initEpidemic(t_Epidemic *pEpidemic, t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel, t_HostStatus *aHostStatus, double *pTotalRate, int epiID)
{
	int retVal, i, t, numToDo, validHosts, thisHost, j;
	int *aHosts;

	fprintf(stdout, "Initialising epidemic %d\n", epiID);
	memset(pEpidemic, 0, sizeof(t_Epidemic));
	*pTotalRate = 0.0;
	retVal = 1;
	/* initialise all host status */
	for (i = 0; i < pHosts->nHosts; i++)
	{
		aHostStatus[i].nGen = _NOT_SET;
		aHostStatus[i].dRate = 0.0;
		aHostStatus[i].eStatus = SUSCEPTIBLE;
	}
	/* do initial infections */
	t = TYPE_I;
	while (retVal && t <= TYPE_II)
	{
		numToDo = pParams->nInitOne;
		if (t == TYPE_II)
		{
			numToDo = pParams->nInitTwo;
		}
		if(numToDo > 0)
		{
			validHosts = 0;
			if (aHosts = malloc(sizeof(int)*pHosts->nHosts))
			{
				for (i = 0; i < pHosts->nHosts; i++)
				{
					if (pHosts->aHosts[i].eType == t)
					{
						aHosts[validHosts] = i;
						validHosts++;
					}
				}
				if (validHosts >= numToDo)
				{
					i = 0;
					while(retVal && i < numToDo)
					{
						thisHost = (int) floor(validHosts*uniformRandom());
						retVal = infectHost(aHosts[thisHost],0.0,_NOT_SET,pEpidemic, pParams, pHosts, pKernel, aHostStatus, pTotalRate);
						if (retVal)
						{
							/* shift the other hosts down one in place to stop a single host being picked twice */
							validHosts--;
							for (j = thisHost; j < validHosts; j++)
							{
								aHosts[j] = aHosts[j + 1]; /* won't run off the end since have already decremented validHosts */
							}
						}
						i++;
					}
				}
				else
				{
					retVal = 0;
				}
			}
			else
			{
				retVal = 0;
			}
		}
		t++;
	}
	return retVal;
}

/*
	debugging function: check all rates are correct given the state
*/
void	checkRates(t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel, t_HostStatus *aHostStatus)
{
	int			i, j;
	double		cacheRate, recalcRate, thisTheta, thisRho;

	for (i = 0; i < pHosts->nHosts; i++)
	{
		cacheRate = aHostStatus[i].dRate;
		recalcRate = _NOT_SET;
		if (aHostStatus[i].eStatus == INFECTED)
		{
			if (pHosts->aHosts[i].eType == TYPE_I)
			{
				recalcRate = pParams->dMuOne;
			}
			else
			{
				recalcRate = pParams->dMuTwo;
			}
		}
		else
		{
			recalcRate = 0.0;
			thisRho = pParams->dRhoOne;
			if (pHosts->aHosts[i].eType == TYPE_II)
			{
				thisRho = pParams->dRhoTwo;
			}
			for (j = 0; j < pHosts->nHosts; j++)
			{
				if (aHostStatus[j].eStatus == INFECTED)
				{
					thisTheta = pParams->dThetaOne;
					if (pHosts->aHosts[j].eType == TYPE_II)
					{
						thisTheta = pParams->dThetaTwo;
					}
					if (aHostStatus[j].nGen >= pParams->nMaxGen)
					{
						thisTheta = 0.0;	/* artificially stop infections once too many generations have passed */
					}
					recalcRate += thisTheta * thisRho * getKernel(i, j, pKernel, pHosts, pParams);
				}
			}
		}
		if (fabs(cacheRate - recalcRate) > 1e-5)
		{
			fprintf(stdout, "%d %f %f\n", i, cacheRate, recalcRate);
		}
	}
}

void dumpEpidemic(t_Params *pParams, t_Hosts *pHosts, t_Epidemic *pEpidemic, FILE *fOut, int itNum, double maxTime)
{
	int  i,j;
	char sTmp[_MAX_STR_LEN];

	/* information on a generation by generation basis */
	if (pParams->eDumpType == DUMP_GENS)
	{
		int *aTypeOneByGen, *aTypeTwoByGen;
		int g;

		aTypeOneByGen = malloc(sizeof(int)*(pParams->nMaxGen + 1));
		if (aTypeOneByGen)
		{
			aTypeTwoByGen = malloc(sizeof(int)*(pParams->nMaxGen + 1));
			if (aTypeTwoByGen)
			{
#ifdef _ONE_LINE_GEN_OUT
				if (itNum == 0)
				{
					fprintf(fOut, "<it>");
				}
				for (g = 0; g <= pParams->nMaxGen; g++)
				{
					if (itNum == 0)
					{
						sprintf(sTmp, ",I_1(%d),I_2(%d)", g, g);
						fprintf(fOut, sTmp);
					}
					if (g > 0)
					{
						fprintf(stdout, "\t");
					}
					sprintf(sTmp, "I_1(%d)\tI_2(%d)", g, g);
					fprintf(stdout, sTmp);
				}
				if (itNum == 0)
				{
					fprintf(fOut, "\n");
				}
				fprintf(stdout, "\n");
				fprintf(fOut, "%d", itNum);
#else
				if (itNum == 0)
				{
					fprintf(fOut, "<it>,<gen>,<n1>,<n2>,<n1+n2>\n");
				}
				fprintf(stdout, "<gen>\t<n1>\t<n2>\t<n1+n2>\n");
#endif
				for (g = 0; g <= pParams->nMaxGen; g++)
				{
					aTypeOneByGen[g] = aTypeTwoByGen[g] = 0;
					for (j = 0; j < pEpidemic->nEntries; j++)
					{
						if (pEpidemic->aEntries[j].nGen == g)
						{
							if (pEpidemic->aEntries[j].eType == TYPE_I)
							{
								aTypeOneByGen[g]++;
							}
							else
							{
								aTypeTwoByGen[g]++;
							}
						}
					}
#ifdef _ONE_LINE_GEN_OUT
					if (g)
					{
						fprintf(stdout, "\t");
					}
					fprintf(stdout, "%d\t%d", aTypeOneByGen[g], aTypeTwoByGen[g]);
					fprintf(fOut, ",%d,%d", aTypeOneByGen[g], aTypeTwoByGen[g]);
#else
					fprintf(stdout, "%d\t%d\t%d\t%d\n", g, aTypeOneByGen[g], aTypeTwoByGen[g], aTypeOneByGen[g] + aTypeTwoByGen[g]);
					fprintf(fOut, "%d,%d,%d,%d,%d\n", itNum, g, aTypeOneByGen[g], aTypeTwoByGen[g], aTypeOneByGen[g] + aTypeTwoByGen[g]);
#endif
				}
#ifdef _ONE_LINE_GEN_OUT
				fprintf(stdout, "\n");
				fprintf(fOut, "\n");
#endif
				free(aTypeTwoByGen);
			}
			free(aTypeOneByGen);
		}
	}
	/* information on a generation by generation basis */
	if (pParams->eDumpType == DUMP_TIMES && maxTime > 0.0)
	{
		double	dStep,thisTime;
		double	aT[N_DUMP_STEPS + 1];
		int		aTypeOneInf[N_DUMP_STEPS + 1];
		int		aTypeTwoInf[N_DUMP_STEPS + 1];

		if (itNum == 0)
		{
			fprintf(fOut, "<it>,<dT>,<n1>,<n2>,<n1+n2>\n");
		}
		fprintf(stdout, "<dT>\t<n1>\t<n2>\t<n1+n2>\n");
		dStep = maxTime / N_DUMP_STEPS;
		for (j = 0; j <= N_DUMP_STEPS; j++)
		{
			thisTime = j*dStep;
			aT[j] = thisTime;
			aTypeOneInf[j] = 0;
			aTypeTwoInf[j] = 0;
			for (i = 0; i < pEpidemic->nEntries; i++)
			{
				if (pEpidemic->aEntries[i].dInfectTime <= thisTime
					&& (pEpidemic->aEntries[i].dRemovalTime < 0 || pEpidemic->aEntries[i].dRemovalTime > thisTime))
				{
					if (pEpidemic->aEntries[i].eType == TYPE_I)
					{
						aTypeOneInf[j]++;
					}
					else
					{
						aTypeTwoInf[j]++;
					}
				}
			}
			fprintf(stdout, "%f\t%d\t%d\t%d\n", aT[j], aTypeOneInf[j], aTypeTwoInf[j], aTypeOneInf[j]+aTypeTwoInf[j]);
			fprintf(fOut, "%d,%f,%d,%d,%d\n", itNum, aT[j], aTypeOneInf[j], aTypeTwoInf[j], aTypeOneInf[j] + aTypeTwoInf[j]);
		}
	}
	if (pParams->bDumpHostStatus)
	{
		char	sDummy[_MAX_STR_LEN],sOutFile[_MAX_STR_LEN];
		char	*pPtr;
		FILE	*fIt;
		int		bInf,nGen;
		double	tI, tR;

		strcpy(sDummy, pParams->sOutFile);
		pPtr = strrchr(sDummy, '.');
		if (pPtr)
		{
			*pPtr = '\0';
		}
		sprintf(sOutFile, "%s_it=%d.csv", sDummy, itNum);
		fIt = fopen(sOutFile, "wb");
		if (fIt)
		{
			fprintf(fIt, "hostID,hostX,hostY,hostType,tI,tR,gen\n");
			for (i = 0; i < pHosts->nHosts; i++)
			{
				bInf = 0;
				j = 0;
				while (bInf == 0 && (j < pEpidemic->nEntries))
				{
					if (pEpidemic->aEntries[j].nHostID == i)
					{
						bInf = 1;
						nGen = pEpidemic->aEntries[j].nGen;
						tI = pEpidemic->aEntries[j].dInfectTime;
						tR = pEpidemic->aEntries[j].dRemovalTime;
					}
					j++;
				}
				if (bInf)
				{
					fprintf(fIt, "%d,%.4f,%.4f,%d,%.4f,%.4f,%d\n", i, pHosts->aHosts[i].dX, pHosts->aHosts[i].dY, pHosts->aHosts[i].eType, tI, tR, nGen);
				}
				else
				{
					fprintf(fIt, "%d,%.4f,%.4f,%d,NA,NA,NA\n", i, pHosts->aHosts[i].dX, pHosts->aHosts[i].dY, pHosts->aHosts[i].eType);
				}
			}
			fclose(fIt);
		}
	}
}

/*
	actually run the epidemics
*/
int runEpidemics(t_Params *pParams, t_Hosts *pHosts, t_Kernel *pKernel)
{
	FILE			*fOut;
	int				*aInfectiveID;
	double			*aInfectiveRate;
	int				retVal, i, j, eventHost, infectingHost, numInfectives, nSteps;
	double			thisExtra, thisTheta, thisRho, runningSum, randDbl, totalRate, timeNow, timeOffset, totalInfectiveRate;
	t_HostStatus	*hostStatus;
	t_Epidemic		sEpidemic;

	retVal = 0;
	fOut = fopen(pParams->sOutFile, "wb");
	if (fOut)
	{
		hostStatus = malloc(sizeof(t_HostStatus) * pHosts->nHosts);
		if (hostStatus)
		{
			aInfectiveID = malloc(sizeof(int) * pHosts->nHosts);
			if (aInfectiveID)
			{
				aInfectiveRate = malloc(sizeof(double) * pHosts->nHosts);
				if (aInfectiveRate)
				{
					retVal = 1;
					i = 0;
					while (retVal && i < pParams->nNumIts)
					{
						/* initialise epidemic */
						timeNow = 0.0;
						retVal = initEpidemic(&sEpidemic, pParams, pHosts, pKernel, hostStatus, &totalRate,i);
						/* run epidemic */
						nSteps = 0;
						while (retVal
								&& (totalRate > 0.0)
								&& (pParams->dMaxTime < 0 || timeNow <= pParams->dMaxTime))
						{
#if 0
							/* check rates every 50 steps (used in debugging) */
							if (nSteps && (nSteps % 50 == 0))
							{
								checkRates(pParams, pHosts, pKernel, hostStatus);
							}
#endif
							/* find time of next event and update current time*/
							randDbl = uniformRandom();
							while (randDbl <= 0.0)
							{
								randDbl = uniformRandom();
							}
							timeOffset = -log(randDbl) / totalRate;
							timeNow = timeNow + timeOffset;

							/* find host that is affected by the event */
							randDbl = totalRate * uniformRandom();
							runningSum = 0.0;
							eventHost = 0;
							do
							{
								runningSum += hostStatus[eventHost].dRate;
								eventHost++;
							} while ((runningSum <= randDbl) && (eventHost < pHosts->nHosts));
							eventHost--;

							/* what happens now depends on whether it is an infection or a recovery */
							if (hostStatus[eventHost].eStatus == SUSCEPTIBLE)
							{
								/* to keep track of generations, need to find which host infected the newly infected one */
								thisRho = pParams->dRhoOne;
								if (pHosts->aHosts[eventHost].eType == TYPE_II)
								{
									thisRho = pParams->dRhoTwo;
								}
								totalInfectiveRate = 0.0;
								numInfectives = 0;
								for (j = 0; j < pHosts->nHosts; j++)
								{
									if (hostStatus[j].eStatus == INFECTED)
									{
										aInfectiveID[numInfectives] = j;
										thisTheta = pParams->dThetaOne;
										if (pHosts->aHosts[j].eType == TYPE_II)
										{
											thisTheta = pParams->dThetaTwo;
										}
										if (hostStatus[j].nGen >= pParams->nMaxGen)
										{
											thisTheta = 0.0;	/* artificially stop infections once too many generations have passed */
										}
										thisExtra = thisTheta * thisRho * getKernel(j, eventHost, pKernel, pHosts, pParams);
										aInfectiveRate[numInfectives] = thisExtra;
										totalInfectiveRate += thisExtra;
										numInfectives++;
									}
								}
								/* find which infected host caused this infection */
								randDbl = totalInfectiveRate * uniformRandom();
								runningSum = 0.0;
								infectingHost = 0;
								do
								{
									runningSum += aInfectiveRate[infectingHost];
									infectingHost++;
								} while ((runningSum <= randDbl) && (infectingHost < numInfectives));
								infectingHost--;
								retVal = infectHost(eventHost, timeNow, aInfectiveID[infectingHost], &sEpidemic, pParams, pHosts, pKernel, hostStatus, &totalRate);
							}
							else
							{
								if (hostStatus[eventHost].eStatus == INFECTED)
								{
									retVal = recoverHost(eventHost, timeNow, &sEpidemic, pParams, pHosts, pKernel, hostStatus, &totalRate);
								}
								else
								{
									fprintf(stderr, "Event triggered by removed host: error\n");
								}
							}
							/*
								since everything has to recover, which puts quite a high lower bound
								on the set of admissible rates, this is a signal of numerical error
								and means the totalRate is really just zero
							*/
							if (totalRate <= _NUMERIC_UNDERFLOW)
							{
								totalRate = 0.0;
							}
#if 0
							{
								double d;

								d = 0;
								for (j = 0; j < pHosts->numHosts; j++)
								{
									d += hostStatus[j].rate;
								}
								fprintf(stdout, "*** %f %f\n", totalRate, d);
							}
#endif
							nSteps++;
						}
						dumpEpidemic(pParams, pHosts, &sEpidemic, fOut, i, pParams->dMaxTime);
						/* add one to iteration number */
						i++;
					}
					free(aInfectiveRate);
				}
				free(aInfectiveID);
			}
			free(hostStatus);
		}
		fclose(fOut);
	}
	return retVal;
}

void freeMemory(t_Hosts *pHosts, t_Kernel *pKernel)
{
	if(pKernel->aKernel)
	{
		free(pKernel->aKernel);
	}
	if (pHosts->nAlloc && pHosts->aHosts)
	{
		free(pHosts->aHosts);
	}
}

/*
	main routine
*/
int main(int argc, char **argv)
{
	t_Params	sParams;
	t_Hosts		sHosts;
	t_Kernel	sKernel;
	int			retVal;

	memset(&sParams, 0, sizeof(t_Params));
	memset(&sHosts, 0, sizeof(t_Hosts));
	memset(&sKernel, 0, sizeof(t_Kernel));
	seedRandom(0); /* passing zero means it uses combination of time and procID as a seed */
#if 0
	if (!setParams(&sParams))
	{
		fprintf(stderr, "Error in setParams()\nExiting\n");
	}
#endif
	if (!(retVal = readParams(&sParams, argc, argv)))
	{
		fprintf(stderr, "Error in readParams()\nExiting\n");
	}
	if (retVal && !(retVal = loadHosts(&sParams, &sHosts)))
	{
		fprintf(stderr, "Error in loadHosts()\nExiting\n");
	}
	if (retVal && !(retVal = calcKernel(&sParams, &sHosts, &sKernel)))
	{
		fprintf(stderr, "Error in calcKernel()\nExiting\n");
	}
	if (retVal && !(retVal = runEpidemics(&sParams, &sHosts, &sKernel)))
	{
		fprintf(stderr, "Error in runEpidemics()\nExiting\n");
	}
	/* free all memory from dynamically allocated structures */
	freeMemory(&sHosts, &sKernel);

	if (retVal == 0)
		return(EXIT_FAILURE);
	return(EXIT_SUCCESS);
}