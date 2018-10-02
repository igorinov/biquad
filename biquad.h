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
double biquad_update(struct iir_filter *filter, double x);
void biquad_init_lowpass(struct iir_filter *filter, double fs, double f);
void biquad_init_bandpass(struct iir_filter *filter, double fs, double f1, double f2);
void biquad_init_bandstop(struct iir_filter *filter, double fs, double f1, double f2);

#endif
