/* Standalone format verifier: loads a group file with the REAL load_group()
 * (same parser dym-mod-metro uses) and checks range + associativity. */
#include "group.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
	if (argc < 2) { fprintf(stderr, "usage: %s groupfile\n", argv[0]); return 2; }
	load_group(argv[1]);   /* aborts on missing identity/inverses */

	/* All table entries must be valid element indices. */
	for (unsigned n = 0; n < P; n++)
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] < 0 || mult[n][m] >= (int)P) {
				printf("FAIL  %-44s mult[%u][%u]=%d out of [0,%u)\n",
				       argv[1], n, m, mult[n][m], P);
				return 1;
			}

	/* Associativity: exhaustive for small groups, sampled for large ones. */
	unsigned long checks = ((unsigned long)P*P*P <= 2000000UL) ? 0 : 300000UL;
	if (checks == 0) {
		for (unsigned a = 0; a < P; a++)
		for (unsigned b = 0; b < P; b++)
		for (unsigned c = 0; c < P; c++)
			if (mult[mult[a][b]][c] != mult[a][mult[b][c]]) {
				printf("FAIL  %-44s non-associative (%u,%u,%u)\n", argv[1], a, b, c);
				return 1;
			}
	} else {
		unsigned long seed = 88172645463325252UL;
		for (unsigned long t = 0; t < checks; t++) {
			seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
			unsigned a = (seed >> 3) % P;
			seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
			unsigned b = (seed >> 3) % P;
			seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
			unsigned c = (seed >> 3) % P;
			if (mult[mult[a][b]][c] != mult[a][mult[b][c]]) {
				printf("FAIL  %-44s non-associative (%u,%u,%u)\n", argv[1], a, b, c);
				return 1;
			}
		}
	}
	printf("PASS  %-44s order=%-5u id=%-4d (id+inv+range+assoc ok)\n", argv[1], P, id);
	return 0;
}
