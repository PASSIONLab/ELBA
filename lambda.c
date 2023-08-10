#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

double alpha = 1.;
double beta = 1.;

int usage(const char *prg)
{
    fprintf(stderr, "Usage: %s [options]\n", prg);
    fprintf(stderr, "Options: -A FLOAT  matching score [%.2f]\n", alpha);
    fprintf(stderr, "         -B FLOAT  mismatch penalty [%.2f]\n", beta);
    fprintf(stderr, "         -h        help message\n");
    return -1;
}

int main(int argc, char *argv[])
{
    int c;

    while ((c = getopt(argc, argv, "A:B:h")) >= 0)
    {
        if      (c == 'A') alpha = atof(optarg);
        else if (c == 'B') beta = atof(optarg);
        else if (c == 'h') return usage(argv[0]);
    }

    beta *= -1;

    double lambda, newlambda, fx, dfx;
    int iter;

    lambda = 0.1;
    for (iter = 0; iter < 100; ++iter)
    {
        if (0.25 * 0.25 * (12 * exp(lambda * beta) + 4 * exp(lambda * alpha)) > 1.) break;
        lambda *= 2.;
    }

    assert(iter != 100);

    for (iter = 0; iter < 100; ++iter)
    {
        fx = 0.25 * 0.25 *(12 * exp(lambda*beta) + 4 * exp(lambda*alpha)) - 1.;

        if (fabs(fx) < 1e-6)
            break;

        dfx = (0.25 * 0.25) * (12 * beta * exp(lambda * beta) + 4 * alpha * exp(lambda * alpha));
        newlambda = lambda - (fx / dfx);
        if (newlambda <= 0) newlambda = 0.000001; /* this shouldn't happen */
        lambda = newlambda;
    }

    assert(iter != 100);

    double target_identity = 0.25 * exp(lambda * alpha);

    printf("target identity = %% %.2f\n", target_identity);

    return 0;
}
