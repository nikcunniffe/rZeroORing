#ifndef _MT19937AR_H_
#define _MT19937AR_H_

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);
/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length);
/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);
/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void);
/* generates a random number on [0,1]-real-interval */
double genrand_real1(void);
/* generates a random number on [0,1)-real-interval */
#ifdef _WIN32
__inline double genrand_real2(void);
#else
inline double genrand_real2(void);
#endif
/* generates a random number on (0,1)-real-interval */
double genrand_real3(void);
/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void);

#endif /* _MT19937AR_H_ */
