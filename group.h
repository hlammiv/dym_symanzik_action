#ifndef GROUP_H
#define GROUP_H

#define PMAX 5040

typedef int group_t;

extern unsigned P;
extern double ReTr[PMAX];
extern double ImTr[PMAX];
extern group_t **mult;   /* P x P Cayley table, allocated contiguously at load
                          * time with stride P (cache-friendly); index mult[n][m] */
extern group_t inv[PMAX];
extern group_t id;

void load_group(const char *);

#endif /* ndef GROUP_H */
