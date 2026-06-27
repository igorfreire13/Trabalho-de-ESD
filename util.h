// =============================================================================
//  util.h - Utilidades visuais e de entrada (terminal)
// =============================================================================
#ifndef UTIL_H
#define UTIL_H
#include "comum.h"

void habilitarANSI();
void printLinha(char c = '-', int n = 57);
long long memoriaRSS_KB();
int lerInteiro(const char *prompt);
double lerDouble(const char *prompt);
void printBarra(double valor, double maximo, int largura = 30, const char *cor = GRN);
void printProgresso(unsigned int atual, unsigned int total);
void printHeader();

#endif
