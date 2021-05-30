#
/*
 *    Copyright (C) 2020
 *    Jan van Katwijk (J.vanKatwijk@gmail.com)
 *    Lazy Chair Computing
 *
 *    This file is part of the dab-cmdline program
 *
 *    dab-cmdline is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    dab-cmdline is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with dab-cmdline; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	"mm_malloc.h"
#include	"viterbi-spiral.h"
#include	<cstring>
#ifdef  __MINGW32__
#include	<intrin.h>
#include	<malloc.h>
#include	<windows.h>
#endif

//
//	It took a while to discover that the polynomes we used
//	in our own "straightforward" implementation was bitreversed!!
//	The official one is on top.
#define K 7
#define POLYS { 0155, 0117, 0123, 0155}
//#define	POLYS	{109, 79, 83, 109}
// In the reversed form the polys look:
//#define POLYS { 0133, 0171, 0145, 0133 }
//#define POLYS { 91, 121, 101, 91 }

#define	METRICSHIFT	0
#define	PRECISIONSHIFT	0
#define	RENORMALIZE_THRESHOLD	137

//
/* ADDSHIFT and SUBSHIFT make sure that the thing returned is a byte. */
#if (K-1<8)
#define ADDSHIFT (8-(K-1))
#define SUBSHIFT 0
#elif (K-1>8)
#define ADDSHIFT 0
#define SUBSHIFT ((K-1)-8)
#else
#define ADDSHIFT 0
#define SUBSHIFT 0
#endif


static uint8_t Partab [] = 
{ 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
  0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0};

//
//	One could create the table above, i.e. a 256 entry
//	odd-parity lookup table by the following function
//	It is now precomputed
void	viterbiSpiral::partab_init (void){
int16_t i,cnt,ti;

	for (i = 0; i < 256; i++){
	   cnt = 0;
	   ti = i;
	   while (ti != 0) {
	      if (ti & 1) cnt++;
	      ti >>= 1;
	   }
	   Partab [i] = cnt & 1;
	}
}

int 	viterbiSpiral::parity (int x){
	/* Fold down to one byte */
	x ^= (x >> 16);
	x ^= (x >> 8);
	return Partab [x];
//	return parityb(x);
}

static inline
void	renormalize (COMPUTETYPE* X, COMPUTETYPE threshold){
int32_t	i;

	if (X [0] > threshold){
	   COMPUTETYPE min = X [0];
	   for (i = 0; i < NUMSTATES; i++)
	      if (min > X[i])
	         min = X[i];
	   for (i = 0; i < NUMSTATES; i++)
	      X[i] -= min;
      }
}

template<typename T>
static inline T* alloc_array(std::size_t align, std::size_t count) {
	T *ptr;
#if defined(__MINGW32__) || defined(_MSC_VER)
	ptr = static_cast<T*>(_aligned_malloc(count * sizeof(T), align));
#else
	if (posix_memalign((void**)&ptr, align, count * sizeof(T)))
		ptr = nullptr;
#endif
	if (ptr == nullptr)
		throw std::bad_alloc();
	return ptr;
}

template<typename T>
static inline void free_array(T *ptr) {
#if defined(__MINGW32__) || defined(_MSC_VER)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

//
	viterbiSpiral::viterbiSpiral (int16_t wordlength) {
int polys [RATE] = POLYS;
int16_t	i, state;
size_t size = wordlength + (K - 1);

	frameBits		= wordlength;
//	partab_init	();

	data = alloc_array<std::remove_reference<decltype(*data)>::type>(ALIGN, size / 8 + 1);
	symbols = alloc_array<std::remove_reference<decltype(*symbols)>::type>(ALIGN, RATE * size);
	vp.decisions = alloc_array<std::remove_reference<decltype(*vp.decisions)>::type>(ALIGN, size);

	for (state = 0; state < NUMSTATES / 2; state++) {
	   for (i = 0; i < RATE; i++)
	      Branchtab [i * NUMSTATES / 2 + state] =
	                     (polys[i] < 0) ^
	                        parity((2 * state) & abs (polys[i])) ? 255 : 0;
	}
//
	init_viterbi (&vp, 0);
}


	viterbiSpiral::~viterbiSpiral	(void) {
		free_array(vp.decisions);
		free_array(data);
		free_array(symbols);
}

static int maskTable [] = {128, 64, 32, 16, 8, 4, 2, 1};
static  inline
uint8_t getbit (uint8_t v, int32_t o) {
        return  (v & maskTable [o]) ? 1 : 0;
}

//static
//uint8_t getbit (uint8_t v, int32_t o) {
//uint8_t	mask	= 1 << (7 - o);
//	return  (v & mask) ? 1 : 0;
//}
	
// depends: POLYS, RATE, COMPUTETYPE
// 	encode was only used for testing purposes
//void encode (/*const*/ unsigned char *bytes, COMPUTETYPE *symbols, int nbits) {
//int	i, k;
//int	polys [RATE] = POLYS;
//int	sr = 0;
//
//// FIXME: this is slowish
//// -- remember about the padding!
//	for (i = 0; i < nbits + (K - 1); i++) {
//	   int b = bytes[i/8];
//	   int j = i % 8;
//	   int bit = (b >> (7-j)) & 1;
//
//	   sr = (sr << 1) | bit;
//	   for (k = 0; k < RATE; k++)
//	      *(symbols++) = parity(sr & polys[k]);
//	}
//}

//	Note that our DAB environment maps the softbits to -127 .. 127
//	we have to map that onto 0 .. 255

void	viterbiSpiral::deconvolve	(int16_t *input, uint8_t *output) {
uint32_t	i;

	init_viterbi (&vp, 0);
	for (i = 0; i < (uint16_t)(frameBits + (K - 1)) * RATE; i ++) {
	   int16_t temp = input [i] + 127;
//	   if (temp < 0) temp = 0;
//	   if (temp > 255) temp = 255;
	   symbols [i] = temp;
	}
//	if (!spiral)
//	   update_viterbi_blk_GENERIC (&vp, symbols, frameBits + (K - 1));
//	else
	   update_viterbi_blk_SPIRAL (&vp, symbols, frameBits + (K - 1));

	chainback_viterbi (&vp, data, frameBits, 0);

	for (i = 0; i < (uint16_t)frameBits; i ++)
	   output [i] = getbit (data [i >> 3], i & 07);
}

/* C-language butterfly */
void	viterbiSpiral::BFLY (int i, int s, COMPUTETYPE * syms,
	                   struct v * vp, decision_t d) {
int32_t j, decision0, decision1;
COMPUTETYPE metric,m0,m1,m2,m3;

	metric =0;
	for (j = 0; j < RATE;j++)
	   metric += (Branchtab [i + j * NUMSTATES/2] ^ syms[s*RATE+j]) >>
	                                                     METRICSHIFT ;
	metric = metric >> PRECISIONSHIFT;
	const COMPUTETYPE max =
	        ((RATE * ((256 - 1) >> METRICSHIFT)) >> PRECISIONSHIFT);

	m0 = vp->old_metrics[i] + metric;
	m1 = vp->old_metrics[i + NUMSTATES / 2] + (max - metric);
	m2 = vp->old_metrics[i] + (max - metric);
	m3 = vp->old_metrics[i + NUMSTATES / 2] + metric;

	decision0 = ((int32_t)(m0 - m1)) > 0;
	decision1 = ((int32_t)(m2 - m3)) > 0;

	vp -> new_metrics[2 * i] = decision0 ? m1 : m0;
	vp -> new_metrics[2 * i + 1] =  decision1 ? m3 : m2;

	d[(2 * i) / decision_bits_per_elem] |=
		(decision0 | decision1 << 1) << ((2 * i) & (decision_bits_per_elem - 1));
}

/* Update decoder with a block of demodulated symbols
 * Note that nbits is the number of decoded data bits, not the number
 * of symbols!
 */
void	viterbiSpiral::update_viterbi_blk_GENERIC (struct v *vp,
					            COMPUTETYPE *syms,
	                                            int16_t nbits){
int32_t  s, i;

	for (s = 0; s < nbits; s++)
	   memset (vp -> decisions + s, 0, sizeof (decision_t));

	for (s = 0; s < nbits; s++){
	   decltype(vp -> old_metrics) tmp;
	   for (i = 0; i < NUMSTATES / 2; i++)
	      BFLY (i, s, syms, vp, vp -> decisions[s]);

	   renormalize (vp -> new_metrics, RENORMALIZE_THRESHOLD);
//     Swap pointers to old and new metrics
	   tmp = vp -> old_metrics;
	   vp -> old_metrics = vp -> new_metrics;
	   vp -> new_metrics = tmp;
	}
}

extern "C" {
#if defined(SSE_AVAILABLE)
void FULL_SPIRAL_sse (int,
#elif defined(NEON_AVAILABLE)
void FULL_SPIRAL_neon (int,
#else
void FULL_SPIRAL_no_sse (int,
#endif
	                 COMPUTETYPE *Y,
	                 COMPUTETYPE *X,
	                 COMPUTETYPE *syms,
	                 uint32_t *dec,
	                 COMPUTETYPE *Branchtab);
}

void	viterbiSpiral::update_viterbi_blk_SPIRAL (struct v *vp,
					           COMPUTETYPE *syms,
					           int16_t nbits){
int32_t s;

	for (s = 0; s < nbits; s++)
	   memset (vp -> decisions + s, 0, sizeof(decision_t));

#if defined(SSE_AVAILABLE)
	FULL_SPIRAL_sse (nbits / 2,
#elif defined(NEON_AVAILABLE)
	FULL_SPIRAL_neon (nbits / 2,
#else
	FULL_SPIRAL_no_sse (nbits / 2,
#endif
	                 vp -> new_metrics,
	                 vp -> old_metrics,
	                 syms,
	                 (uint32_t*)vp -> decisions, Branchtab);
}

//
/* Viterbi chainback */
void	viterbiSpiral::chainback_viterbi (struct v *vp,
	                            uint8_t *data, /* Decoded output data */
	                            int16_t nbits, /* Number of data bits */
	                            uint16_t endstate){ /*Terminal encoder state */
/* Make room beyond the end of the encoder register so we can
 * accumulate a full byte of decoded data
 */
	endstate = (endstate % NUMSTATES) << ADDSHIFT;
/* The store into data[] only needs to be done every 8 bits.
 * But this avoids a conditional branch, and the writes will
 * combine in the cache anyway
 */
	while (nbits-- != 0){
	   int k;
//	   int l	= (endstate >> ADDSHIFT) / 32;
//	   int m	= (endstate >> ADDSHIFT) % 32;
	   k = (vp -> decisions [(K - 1) + nbits] [(endstate >> ADDSHIFT) / decision_bits_per_elem] >>
	                       ((endstate>>ADDSHIFT) % decision_bits_per_elem)) & 1;
	   endstate = (endstate >> 1) | (k << (K - 2 + ADDSHIFT));
	   data [nbits >> 3] = endstate >> SUBSHIFT;
	}
}

/* Initialize Viterbi decoder for start of new frame */
void 	viterbiSpiral::init_viterbi (struct v *p, int16_t starting_state){
struct v *vp = p;
int32_t i;

	for (i = 0; i < NUMSTATES; i++)
	   vp -> metrics1[i] = 63;

	vp -> old_metrics = vp -> metrics1;
	vp -> new_metrics = vp -> metrics2;
/* Bias known start state */
	vp -> old_metrics[starting_state & (NUMSTATES-1)] = 0;
}

