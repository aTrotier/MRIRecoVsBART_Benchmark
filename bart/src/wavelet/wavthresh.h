
#include <stdbool.h>

extern const struct operator_p_s* prox_wavelet_thresh_create(unsigned int N, const long dims[N], unsigned int flags, unsigned int jflags, const long minsize[N], float lambda, bool randshift);


extern void wavthresh_rand_state_set(const struct operator_p_s* op, int x);


