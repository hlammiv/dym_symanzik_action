#ifndef GROUP_H
#define GROUP_H

#define PMAX 1080

typedef int group_t;

extern unsigned P;
extern double ReTr[PMAX];
extern double ImTr[PMAX];
extern group_t mult[PMAX][PMAX];
extern group_t inv[PMAX];
extern group_t id;

void load_group(const char *);

#endif /* ndef GROUP_H */
