// =============================================================================
//  comum.h - Definicoes compartilhadas: constantes, estruturas e globais
// =============================================================================
#ifndef COMUM_H
#define COMUM_H

#include <cstdint>

// ---- Tamanhos das estruturas ------------------------------------------------
#define TAMANHO_TABELA 500009
#define TAM_HASH_OTM 600011
#define K 2
#define TAM_BLOOM 5000000
#define NUM_HASHES 3
#define LIMITE_MAX_AVL 100000

// ---- Macros de bit para Bloom Filter ----------------------------------------
#define SET_BIT(arr, k) ((arr)[(k) / 8] |= (1 << ((k) % 8)))
#define TEST_BIT(arr, k) ((arr)[(k) / 8] & (1 << ((k) % 8)))

// ---- Cores ANSI -------------------------------------------------------------
#define RST "\033[0m"
#define BOLD "\033[1m"
#define RED "\033[1;31m"
#define GRN "\033[1;32m"
#define YEL "\033[1;33m"
#define BLU "\033[1;34m"
#define CYN "\033[1;36m"
#define WHT "\033[1;37m"
#define DIM "\033[2m"

// ---- Estruturas de dados ----------------------------------------------------
typedef struct NoLista
{
    double bounding_area;
    unsigned int sample_id;
    struct NoLista *prox;
} NoLista;

// Hash otimizada: Open Addressing + Tombstone
typedef struct
{
    unsigned int sample_id;
    double bounding_area;
    bool ocupado, removido;
} HashSlot;
typedef struct
{
    HashSlot *tabela;
} TabelaHashOtm;

typedef struct NoAVL
{
    double bounding_area;
    unsigned int sample_id;
    char frame_name[256];
    struct NoAVL *esq, *dir;
    int altura;
} NoAVL;

typedef struct NoHash
{
    unsigned int sample_id;
    double bounding_area;
    struct NoHash *prox;
} NoHash;
typedef struct
{
    NoHash *baldes[TAMANHO_TABELA];
    long long total_colisoes;
} TabelaHash;

typedef struct NoKD
{
    double ponto[K];
    unsigned int sample_id;
    int frame_id;
    struct NoKD *esq, *dir;
} NoKD;

typedef struct
{
    uint8_t *bits;
    int tamanho;
} BloomFilter;

// Resumo estatistico de uma amostra
typedef struct
{
    double minimo, q1, mediana, media, q3, maximo, desvio, iqr;
    double iqr_inf, iqr_sup;
    double zsc_inf, zsc_sup;
} ResumoEstatistico;

// Ponto coletado da KD-Tree (x, y, id, frame)
struct PontoKD
{
    double x, y;
    unsigned int id;
    int frame;
};

// ---- Globais de restricao (definidos em main.cpp) ---------------------------
extern int g_total_avl_nodes;
extern bool g_janela_aberta;
extern int g_limite_avl;          // R2: teto de nos da AVL (runtime)
extern long long g_avl_rejeitados;
extern volatile bool g_sink_b;    // sinks: impedem o compilador de eliminar loops
extern volatile double g_sink_d;

// Guard RAII: salva e restaura o estado global da AVL ao sair de escopo
struct EstadoAVLGuard
{
    int total;
    int limite;
    long long rej;
    EstadoAVLGuard() : total(g_total_avl_nodes), limite(g_limite_avl), rej(g_avl_rejeitados) {}
    ~EstadoAVLGuard()
    {
        g_total_avl_nodes = total;
        g_limite_avl = limite;
        g_avl_rejeitados = rej;
    }
};

#endif
