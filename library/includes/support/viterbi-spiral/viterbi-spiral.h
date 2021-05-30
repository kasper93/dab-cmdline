#
#ifndef	__VITERBI_SPIRAL__
#define	__VITERBI_SPIRAL__
/*
 * 	Viterbi.h according to the SPIRAL project
 */
#include	"dab-constants.h"

//	For our particular viterbi decoder, we have
#define	RATE	4
#define NUMSTATES 64
#define ALIGN 16
#define COMPUTETYPE uint32_t

//decision_t is a BIT vector
using decision_t = uint32_t[NUMSTATES / (sizeof(uint32_t) * 8)];
using metric_t = COMPUTETYPE[NUMSTATES];

static constexpr size_t decision_bits_per_elem =
	sizeof(std::remove_all_extents<decision_t>::type) * 8;

/* State info for instance of Viterbi decoder
 */

struct v {
/* path metric buffer 1 */
	alignas(ALIGN) metric_t metrics1;
/* path metric buffer 2 */
	alignas(ALIGN) metric_t metrics2;
/* Pointers to path metrics, swapped on every bit */
	COMPUTETYPE *old_metrics, *new_metrics;
	decision_t *decisions;   /* decisions */
};

class	viterbiSpiral {
public:
		viterbiSpiral	(int16_t);
		~viterbiSpiral	(void);
	void	deconvolve	(int16_t *, uint8_t *);
private:

	struct v	vp;
	alignas(ALIGN) COMPUTETYPE Branchtab	[NUMSTATES / 2 * RATE];
//	int	parityb		(uint8_t);
	int	parity		(int);
	void	partab_init	(void);
//	uint8_t	Partab	[256];
	void	init_viterbi	(struct v *, int16_t);
	void	update_viterbi_blk_GENERIC	(struct v *, COMPUTETYPE *,
	                                         int16_t);
	void	update_viterbi_blk_SPIRAL	(struct v *, COMPUTETYPE *,
	                                         int16_t);
	void	chainback_viterbi (struct v *, uint8_t *, int16_t, uint16_t);
	struct v *viterbi_alloc (int32_t);
	void	BFLY		(int32_t, int, COMPUTETYPE *,
	                         struct v *, decision_t);
//	uint8_t *bits;
	uint8_t *data;
	COMPUTETYPE *symbols;
	int16_t	frameBits;
};

#endif

