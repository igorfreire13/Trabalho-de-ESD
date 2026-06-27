// =============================================================================
//  benchmark.h - Benchmarks, escalabilidade e teste de restricao R2
// =============================================================================
#ifndef BENCHMARK_H
#define BENCHMARK_H
#include "comum.h"

void executarBenchmark(NoLista *lista, TabelaHashOtm *hotm, NoAVL *avl,
                       TabelaHash *hash, BloomFilter *bloom, unsigned int id_ref);
void testeEscalabilidade();
void testarRestricaoR2();

#endif
