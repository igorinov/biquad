#include <math.h>
#include <stdlib.h>

#include "biquad.h"

typedef struct {
    float re;
    float im;
} complex_s;

typedef struct {
    double re;
    double im;
} complex_d;

iir_filter_t *biquad_create(int sections)
{
    iir_filter_t *filter;
    int sect_ord = 2;

    filter = (iir_filter_t *) calloc(1, sizeof(iir_filter_t));
    if (!filter)
        return NULL;

    filter->sections = sections;
    filter->sect_ord = sect_ord;

    filter->a = (double *) calloc(sections * (sect_ord + 1), sizeof(double));
    filter->b = (double *) calloc(sections * (sect_ord + 1), sizeof(double));
    filter->d = (double *) calloc((sections + 1) * sect_ord, sizeof(double));

    if (!filter->a || !filter->b || !filter->d) {
        biquad_delete(filter);
        return NULL;
    }

    return filter;
}

void biquad_delete(iir_filter_t *filter)
{
    if (!filter)
        return;

    if (filter->a) {
        free(filter->a);
        filter->a = NULL;
    }

    if (filter->b) {
        free(filter->b);
        filter->b = NULL;
    }

    if (filter->d) {
        free(filter->d);
        filter->d = NULL;
    }

    free(filter);
}

double biquad_update(struct iir_filter *filter, double x)
{
    double *a = filter->a;
    double *b = filter->b;
    double *d = filter->d;
    int stages = filter->sections;
    int k;
    double y;

    for (k = 0; k < stages; k += 1) {
        y = x * b[0];
        y += d[0] * b[1] + d[1] * b[2];
        y -= d[2] * a[1] + d[3] * a[2];
        d[1] = d[0];
        d[0] = x;
        x = y;
        d += 2;
        a += 3;
        b += 3;
    }
    d[1] = d[0];
    d[0] = x;

    return x;
}

void biquad_init_lowpass(struct iir_filter *filter, double fs, double f) {
    double *a = filter->a;
    double *b = filter->b;
    double w = 2 * M_PI * f / fs;
    double phi, alpha;
    int i, k;
    int n;

    n = filter->sections;
    for (i = 0; i < n; i += 1) {
        k = n - i - 1;
        phi = M_PI / (4.0 * n) * (k * 2 + 1);
        alpha = sin(w) * cos(phi);

        b[0] = (1.0 - cos(w)) / (2 * (1.0 + alpha));
        b[1] = (1 - cos(w)) / (1.0 + alpha);
        b[2] = (1.0 - cos(w)) / (2 * (1.0 + alpha));
        a[0] = 1;
        a[1] = -2 * cos(w) / (1.0 + alpha);
        a[2] = (1.0 - alpha) / (1.0 + alpha);
        a += 3;
        b += 3;
    }

    for (i = 0; i < (n + 1) * 2; i += 1)
        filter->d[i] = 0;
}

void biquad_init_highpass(struct iir_filter *filter, double fs, double f) {
    double *a = filter->a;
    double *b = filter->b;
    double w = 2 * M_PI * f / fs;
    double phi, alpha;
    int n, i, k;

    n = filter->sections;

    for (i = 0; i < n; i += 1) {
        k = n - i - 1;
        phi = M_PI / (4.0 * n) * (k * 2 + 1);
        alpha = sin(w) * cos(phi);

        b[0] = (1.0 + cos(w)) / (2 * (1.0 + alpha));
        b[1] = -(1 + cos(w)) / (1.0 + alpha);
        b[2] = (1.0 + cos(w)) / (2 * (1.0 + alpha));
        a[1] = -2 * cos(w) / (1.0 + alpha);
        a[2] = (1.0 - alpha) / (1.0 + alpha);
        a += 3;
        b += 3;
    }

    for (i = 0; i < (n + 1) * 2; i += 1)
        filter->d[i] = 0;
}

static void complex_square(complex_d *s)
{
    double x, y;

    x = s->re;
    y = s->im;

    s->re = x * x - y * y;
    s->im = 2 * x * y;
}

static void complex_sqrt(complex_d *s)
{
    double x, y;
    double r, phi;

    x = s->re;
    y = s->im;
    if (x == 0 && y == 0) {
        return;
    }

    /* Converting to polar */
    phi = atan2(y, x);
    r = sqrt(x * x + y * y);

    /* Square root */
    phi /= 2;
    r = sqrt(r);

    /* Back to cartesian */
    s->re = r * cos(phi);
    s->im = r * sin(phi);
}

static void complex_mul(complex_d *a, complex_d *b)
{
    double x, y;

    x = a->re * b->re - a->im * b->im;
    y = a->re * b->im + a->im * b->re;

    a->re = x;
    a->im = y;
}

static void complex_div(complex_d *p, complex_d *q)
{
    double x, y;
    double r, phi;
    complex_d b;

    x = q->re;
    y = q->im;
    if (x == 0 && y == 0) {
        return;
    }

    /* Converting to polar */
    phi = atan2(y, x);
    r = sqrt(x * x + y * y);

    /* Reciprocal */
    phi = -phi;
    r = 1.0 / r;

    /* Back to cartesian */
    b.re = r * cos(phi);
    b.im = r * sin(phi);

    complex_mul(p, &b);
}

void biquad_init_bandpass(struct iir_filter *filter, double fs, double f1, double f2)
{
    double ts = 1.0 / fs;
    double bw, f;
    double w = 2 * M_PI * f / fs;
    complex_d p, q;
    double phi;
    complex_d z, p_lp, p_bp;
    double k, x, y;
    double wa1, wa2, wa;
    double *a = filter->a;
    double *b = filter->b;
    int n, i;

    /* Map to continuous-time frequencies (pre-warp) */

    wa1 = 2 * fs * tan(M_PI * f1 * ts);
    wa2 = 2 * fs * tan(M_PI * f2 * ts);

    bw = wa2 - wa1;
    wa = sqrt(wa1 * wa2);

    n = filter->sections;

    for (i = 0; i < n; i += 1) {
        phi = M_PI / 2 + M_PI * (2 * i + 1) / (n * 2);
        x = cos(phi);
        y = sin(phi);

        p_lp.re = x * bw / (wa * 2);
        p_lp.im = y * bw / (wa * 2);

        /* Map every low-pass pole to a pair of band-bass poles */

        p_bp = p_lp;
        complex_square(&p_bp);
        p_bp.re = 1 - p_bp.re;
        p_bp.im = 0 - p_bp.im;
        complex_sqrt(&p_bp);
        x = p_lp.re - p_bp.im;
        y = p_lp.im + p_bp.re;
        x *= M_PI * f * ts;
        y *= M_PI * f * ts;

        /* Bilinear transform */

        p.re = 1.0 + x;
        p.im = y;
        q.re = 1.0 - x;
        q.im = -y;
        complex_div(&p, &q);
        x = p.re;
        y = p.im;

        b[0] = 1;
        b[1] = 0;
        b[2] = -1;
        a[0] = 1;
        a[1] = -2 * x;
        a[2] = x * x + y * y;

        /* Scale the parameters to get unity gain in the bassband */

        z.re = cos(w);
        z.im = sin(w);

        p.re = b[2];
        p.im = 0;
        complex_div(&p, &z);
        p.re += b[1];
        complex_div(&p, &z);
        p.re += b[0];

        q.re = a[2];
        q.im = 0;
        complex_div(&q, &z);
        q.re += a[1];
        complex_div(&q, &z);
        q.re += 1;

        complex_div(&p, &q);

        x = p.re;
        y = p.im;
        k = 1.0 / sqrt(x * x + y * y);

        b[0] *= k;
        b[1] *= k;
        b[2] *= k;

        a += filter->sect_ord + 1;
        b += filter->sect_ord + 1;
    }

    for (i = 0; i < (n + 1) * 2; i += 1)
        filter->d[i] = 0;
}
