// =============================================================================
//  estruturas.h - Estruturas de dados: Lista, AVL, Hash, KD-Tree, Bloom
// =============================================================================
#ifndef ESTRUTURAS_H
#define ESTRUTURAS_H
#include "comum.h"
#include <vector>

// ---- Lista encadeada --------------------------------------------------------
NoLista *inserirLista(NoLista *cab, double area, unsigned int id);
bool buscarLista(NoLista *cab, double area);
NoLista *removerLista(NoLista *cab, unsigned int id);
void liberarLista(NoLista *cab);

// ---- Hash otimizada (open addressing) ---------------------------------------
TabelaHashOtm *criarHashOtm();
void inserirHashOtm(TabelaHashOtm *h, unsigned int id, double area);
bool buscarHashOtm(TabelaHashOtm *h, unsigned int id, double *area);
void removerHashOtm(TabelaHashOtm *h, unsigned int id);
void liberarHashOtm(TabelaHashOtm *h);

// ---- Arvore AVL -------------------------------------------------------------
int alturaAVL(NoAVL *n);
NoAVL *criarNoAVL(double area, unsigned int id, const char *frame);
NoAVL *inserirAVL(NoAVL *n, double area, unsigned int id, const char *frame);
NoAVL *minAVL(NoAVL *n);
NoAVL *removerAVL(NoAVL *r, double area, unsigned int id);
NoAVL *buscarAVL(NoAVL *r, double area);
void liberarAVL(NoAVL *r);
NoAVL *podarAVL(NoAVL *r, double mn, double mx);
void coletarAVL(NoAVL *r, std::vector<double> &vals);

// ---- Hash encadeamento ------------------------------------------------------
int funcaoHash(unsigned int id);
TabelaHash *criarTabelaHash();
void inserirHash(TabelaHash *t, unsigned int id, double area);
bool buscarHash(TabelaHash *t, unsigned int id, double *area);
bool buscarHashRestrita(TabelaHash *t, unsigned int id, double *area);
void removerHash(TabelaHash *t, unsigned int id);
void liberarTabelaHash(TabelaHash *t);

// ---- KD-Tree ----------------------------------------------------------------
NoKD *criarNoKD(double pt[], unsigned int id, int frame);
NoKD *inserirKD(NoKD *r, double pt[], unsigned int id, int frame);
bool buscarKD(NoKD *r, double pt[]);
void liberarKD(NoKD *r);

// ---- Bloom Filter -----------------------------------------------------------
BloomFilter *criarBloom();
void inserirBloom(BloomFilter *f, unsigned int id);
bool buscarBloom(BloomFilter *f, unsigned int id);
void liberarBloom(BloomFilter *f);

#endif
