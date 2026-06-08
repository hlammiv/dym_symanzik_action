#include "group.h"

#include <stdio.h>
#include <stdlib.h>

unsigned P;
double ReTr[PMAX];
double ImTr[PMAX];
group_t mult[PMAX][PMAX];
group_t inv[PMAX];
group_t id;

void load_group(const char *fn) {
	FILE *fin = fopen(fn, "r");
	/* File format for group specification is as follows:
	 *
	 *  <order>
	 *  ReTr0 ReTr1 ReTr2 ...
	 *  ImTr0 ImTr1 ImTr2 ...
	 *  0x0 0x1 0x2 ...
	 *  1x0 1x1 1x2 ...
	 *  ...
	 *
	 * For a group of order P, there are thus 1 + P + P^2 entries. Whitespace is ignored.
	 */
	fscanf(fin, "%d", &P);
	if (P > PMAX) {
		fprintf(stderr, "Order of group too large: %d > %d\n", P, PMAX);
		abort();
	}
	for (unsigned n = 0; n < P; n++)
		fscanf(fin, "%lf", &ReTr[n]);
        for (unsigned n = 0; n < P; n++)
                fscanf(fin, "%lf", &ImTr[n]);
	for (unsigned n = 0; n < P; n++)
		for (unsigned m = 0; m < P; m++)
			fscanf(fin, "%d", &mult[n][m]);
	fclose(fin);

	// Find the identity.
	char id_found = 0;
	for (unsigned n = 0; n < P; n++) {
		char is_id = 1;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] != m || mult[m][n] != m)
				is_id = 0;
		if (is_id) {
			id_found = 1;
			id = n;
			break;
		}
	}
	if (!id_found) {
		fprintf(stderr, "Group does not have identity element\n");
		abort();
	}

	// Find inverses.
	for (unsigned n = 0; n < P; n++) {
		char inv_found = 0;
		for (unsigned m = 0; m < P; m++)
			if (mult[n][m] == id && mult[m][n] == id) {
				inv[n] = m;
				inv_found = 1;
				break;
			}
		if (!inv_found) {
			fprintf(stderr, "Group does not have inverses\n");
			abort();
		}
	}
}
