#ifndef __biquad_h
#define __biquad_h

typedef struct iir_filter {
    int sections;
    int sect_ord;
    double *a;
    double *b;
    double *d;
} iir_filter_t;

iir_filter_t *biquad_create(int sections);
void biquad_delete(iir_filter_t *filter);
double biquad_update(iir_filter_t *filter, double x);
void iir_freq_resp(iir_filter_t *filter, double *h, double fs, double f);
void biquad_zero(struct iir_filter *filter);
void biquad_init_lowpass(iir_filter_t *filter, double fs, double f);
void biquad_init_highpass(iir_filter_t *filter, double fs, double f);
void biquad_init_bandpass(iir_filter_t *filter, double fs, double f1, double f2);
void biquad_init_bandstop(iir_filter_t *filter, double fs, double f1, double f2);

#endif
