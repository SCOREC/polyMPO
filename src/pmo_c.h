#ifndef PMO_C_H
extern "C" {
void polympo_initialize();
void polympo_finalize();
void polympo_modifyArray(int size, double* array);
}
#endif
