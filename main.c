#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// Convert a string number in given base to decimal
long long toDecimal(const char *num, int base) {
    long long result = 0;
    for (int i = 0; num[i]; i++) {
        int digit;
        if (isdigit(num[i])) digit = num[i] - '0';
        else if (isalpha(num[i])) digit = tolower(num[i]) - 'a' + 10;
        else continue;
        result = result * base + digit;
    }
    return result;
}

// Solve k x k system using Gaussian elimination
void solvePolynomialCoefficients(double **mat, double *vals, int k, double *coeff) {
    int i, j, l;
    double factor;

    // Gaussian elimination
    for (i = 0; i < k; i++) {
        // Pivot
        for (l = i + 1; l < k; l++) {
            if (fabs(mat[i][i]) < fabs(mat[l][i])) {
                for (j = 0; j < k; j++) {
                    double tmp = mat[i][j]; mat[i][j] = mat[l][j]; mat[l][j] = tmp;
                }
                double tmpb = vals[i]; vals[i] = vals[l]; vals[l] = tmpb;
            }
        }
        // Eliminate
        for (j = i + 1; j < k; j++) {
            factor = mat[j][i] / mat[i][i];
            for (int m = i; m < k; m++) mat[j][m] -= factor * mat[i][m];
            vals[j] -= factor * vals[i];
        }
    }

    // Back-substitution
    for (i = k - 1; i >= 0; i--) {
        coeff[i] = vals[i];
        for (j = i + 1; j < k; j++) coeff[i] -= mat[i][j] * coeff[j];
        coeff[i] /= mat[i][i];
    }
}

int main() {
    // Example input
    int n = 4; // number of roots
    int k = 3; // minimum roots required (degree m = k-1)
    int bases[] = {10, 2, 10, 4};
    const char *values[] = {"4", "111", "12", "213"};

    // Convert roots to decimal
    long long decimalRoots[n];
    printf("Decimal Roots:\n");
    for (int i = 0; i < n; i++) {
        decimalRoots[i] = toDecimal(values[i], bases[i]);
        printf("Root %d: %lld\n", i + 1, decimalRoots[i]);
    }

    // Build Vandermonde matrix for first k roots
    double **mat = (double **)malloc(k * sizeof(double *));
    for (int i = 0; i < k; i++) {
        mat[i] = (double *)malloc(k * sizeof(double));
        mat[i][0] = 1;
        for (int j = 1; j < k; j++) mat[i][j] = mat[i][j - 1] * decimalRoots[i];
    }

    // Polynomial values at roots = 0
    double vals[k];
    for (int i = 0; i < k; i++) vals[i] = 0;

    // Solve for coefficients
    double coeff[k];
    solvePolynomialCoefficients(mat, vals, k, coeff);

    printf("\nPolynomial Coefficients:\n");
    for (int i = 0; i < k; i++) printf("a%d = %.0lf\n", i, coeff[i]);

    // Free memory
    for (int i = 0; i < k; i++) free(mat[i]);
    free(mat);

    return 0;
}