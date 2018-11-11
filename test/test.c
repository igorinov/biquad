#include <math.h>
#include <stdio.h>
#include <string.h>

#include "biquad.h"

static inline double _power(double x, double y)
{
    return x * x + y * y;
}

double eps = 0.000000001;

/* maximum allowed filter gain: 1 + eps1 */
double eps1 = 0.000000000001;

int check_response(iir_filter_t *filter, double fs, double *ft, double *p, int n)
{
    double h[2];
    double f, x, df;
    int i;
    int errors = 0;
    int m = 4096;

    for (i = 0; i < n; i += 1) {
        f = ft[i];
        iir_freq_resp(filter, h, fs, f);
        x = _power(h[0], h[1]);
        if (fabs(x - p[i]) > eps) {
            printf("    H(%.3lf) = %.9lf, expected %.9lf\n", f, x, p[i]);
            errors += 1;
        }
    }

    df = fs / (2.0 * m);
    for (i = 0; i < m; i += 1) {
        f = i * df;
        iir_freq_resp(filter, h, fs, f);
        x = _power(h[0], h[1]);
        if (x - 1 > eps1) {
            printf("    H(%.3lf) = %.9lf > 1\n", f, x);
            errors += 1;
            break;
        }
    }

    return errors;
}

int main(int argc, char **argv)
{
    iir_filter_t *filter;
    double fs = 8000;
    double f0 = 880;
    double fl, fu;
    double k;
    double test_f[5];
    double test_p[5];
    int errors = 0;
    int n;

    for (n = 2; n <= 7; n += 1) {
        printf("Sections: %d\n", n);
        filter = biquad_create(n);

        test_f[0] = 0;
        test_f[1] = f0;
        test_f[2] = fs / 2;

        printf("Low-pass filter frequency response test\n");
        biquad_init_lowpass(filter, fs, f0);

        test_p[0] = 1.0;
        test_p[1] = 0.5;
        test_p[2] = 0;

        errors += check_response(filter, fs, test_f, test_p, 3);

        printf("Highpass filter frequency response test\n");
        biquad_init_highpass(filter, fs, f0);

        test_p[0] = 0;
        test_p[1] = 0.5;
        test_p[2] = 1.0;

        errors += check_response(filter, fs, test_f, test_p, 3);

        fl = f0 * pow(2, -1 / 12.0);
        fu = f0 * pow(2, 1 / 12.0);

        test_f[0] = 0;
        test_f[1] = fl;
        test_f[2] = f0;
        test_f[3] = fu;
        test_f[4] = fs / 2;

        printf("Bandpass filter frequency response test\n");
        biquad_init_bandpass(filter, fs, fl, fu);

        test_p[0] = 0;
        test_p[1] = 0.5;
        test_p[2] = 1.0;
        test_p[3] = 0.5;
        test_p[4] = 0;

        errors += check_response(filter, fs, test_f, test_p, 5);

        printf("Bandstop filter frequency response test\n");
        biquad_init_bandstop(filter, fs, fl, fu);

        test_p[0] = 1.0;
        test_p[1] = 0.5;
        test_p[2] = 0;
        test_p[3] = 0.5;
        test_p[4] = 1.0;

        errors += check_response(filter, fs, test_f, test_p, 5);

        biquad_delete(filter);
    }

    if (errors == 0)
        printf("All tests passed.\n");
    else
        printf("Failed tests: %d\n", errors);

    return 0;
}

