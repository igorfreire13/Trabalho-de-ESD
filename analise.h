// =============================================================================
//  analise.h - Operacoes adicionais: estatistica, histograma, analise
//  espacial (KD-Tree) e ingestao do dataset
// =============================================================================
#ifndef ANALISE_H
#define ANALISE_H
#include "comum.h"
#include <vector>
#include <string>

// ---- Estatistica ------------------------------------------------------------
double calcularMedia(NoLista *cab, int *total);
double calcularDP(NoLista *cab, double media, int total);
ResumoEstatistico calcularResumo(std::vector<double> vals);

// ---- Histograma -------------------------------------------------------------
void exibirHistograma(NoLista *cab, NoAVL *avl);

// ---- Analise espacial (KD-Tree) ---------------------------------------------
void analiseEspacialKD(NoKD *kd);

// ---- Dataset ----------------------------------------------------------------
void gerarDataset(const std::string &nome, int n);
void carregarDataset(NoLista **lista, TabelaHashOtm *hotm, NoAVL **avl,
                     TabelaHash *hash, NoKD **kd, BloomFilter *bloom);

#endif
