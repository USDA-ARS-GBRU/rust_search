#include <stdio.h>
#include <stddef.h>
#include "primer3/src/thal.h"

int main() {
    printf("%zu\n", offsetof(thal_results, temp));
    return 0;
}
