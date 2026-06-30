// =============================================================================
//  SISTEMA DE ANALISE DE DETECCOES DE UVA - arquivo unico (main.cpp)
//  Estruturas de Dados aplicadas a Visao Computacional
// =============================================================================


// ----- (origem: comum.h) -----

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


// ----- (origem: util.h) -----

// =============================================================================
//  util.h - Utilidades visuais e de entrada (terminal)
// =============================================================================
#ifndef UTIL_H
#define UTIL_H

void habilitarANSI();
void printLinha(char c = '-', int n = 57);
long long memoriaRSS_KB();
int lerInteiro(const char *prompt);
double lerDouble(const char *prompt);
void printBarra(double valor, double maximo, int largura = 30, const char *cor = GRN);
void printProgresso(unsigned int atual, unsigned int total);
void printHeader();

#endif


// ----- (origem: estruturas.h) -----

// =============================================================================
//  estruturas.h - Estruturas de dados: Lista, AVL, Hash, KD-Tree, Bloom
// =============================================================================
#ifndef ESTRUTURAS_H
#define ESTRUTURAS_H
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


// ----- (origem: analise.h) -----

// =============================================================================
//  analise.h - Operacoes adicionais: estatistica, histograma, analise
//  espacial (KD-Tree) e ingestao do dataset
// =============================================================================
#ifndef ANALISE_H
#define ANALISE_H
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


// ----- (origem: benchmark.h) -----

// =============================================================================
//  benchmark.h - Benchmarks, escalabilidade e teste de restricao R2
// =============================================================================
#ifndef BENCHMARK_H
#define BENCHMARK_H

void executarBenchmark(NoLista *lista, TabelaHashOtm *hotm, NoAVL *avl,
                       TabelaHash *hash, BloomFilter *bloom, unsigned int id_ref);
void testeEscalabilidade();
void testarRestricaoR2();

#endif


// ----- (origem: util.cpp) -----

#include <cstdio>
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
extern "C" BOOL WINAPI K32GetProcessMemoryInfo(HANDLE, PROCESS_MEMORY_COUNTERS *, DWORD);
#endif
using namespace std;

void habilitarANSI()
{
#ifdef _WIN32
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD m = 0;
    if (GetConsoleMode(h, &m))
        SetConsoleMode(h, m | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
    SetConsoleOutputCP(CP_UTF8);
#endif
}

void printLinha(char c, int n)
{
    printf(DIM);
    for (int i = 0; i < n; i++)
        putchar(c);
    printf(RST "\n");
}

// Retorna a memoria RAM (Working Set / RSS) realmente usada pelo processo, em KB.
long long memoriaRSS_KB()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return (long long)(pmc.WorkingSetSize / 1024);
#endif
    return -1;
}

// --- Leitura robusta de input: nunca trava o menu com entrada invalida ---
int lerInteiro(const char *prompt)
{
    int v;
    while (true)
    {
        printf("%s", prompt);
        if (cin >> v)
            return v;
        cin.clear();
        cin.ignore(10000, '\n');
        printf(RED "  Entrada invalida. Digite um numero inteiro.\n" RST);
    }
}

double lerDouble(const char *prompt)
{
    double v;
    while (true)
    {
        printf("%s", prompt);
        if (cin >> v)
            return v;
        cin.clear();
        cin.ignore(10000, '\n');
        printf(RED "  Entrada invalida. Digite um numero.\n" RST);
    }
}

// Barra de progresso horizontal  [####....] pct%
void printBarra(double valor, double maximo, int largura, const char *cor)
{
    int cheio = (maximo > 0) ? (int)(valor / maximo * largura) : 0;
    if (cheio > largura)
        cheio = largura;
    printf("%s[", cor);
    for (int i = 0; i < largura; i++)
        putchar(i < cheio ? '#' : '.');
    printf("]" RST);
}

// Progresso de ingestão -reescreve a linha com \r
void printProgresso(unsigned int atual, unsigned int total)
{
    int pct = (total > 0) ? (int)((double)atual / total * 100) : 0;
    int blocos = pct * 40 / 100;
    printf("\r  " GRN "[");
    for (int i = 0; i < 40; i++)
        putchar(i < blocos ? '#' : '.');
    printf("]" RST " %3d%%  (%u amostras)", pct, atual);
    fflush(stdout);
}

// Cabeçalho do sistema
void printHeader()
{
    printLinha('=');
    printf(CYN BOLD
           "   SISTEMA DE ANALISE DE DETECCOES DE UVA\n"
           "   Estruturas de Dados aplicadas a Visao Computacional - 2026\n" RST);
    printLinha('=');
}


// ----- (origem: estruturas.cpp) -----

#include <cstdio>
#include <cstring>
#include <cmath>
using namespace std;

// =============================================================================
//  2. LISTA ENCADEADA
// -----------------------------------------------------------------------------
//  Baseline O(n): estrutura mais simples (nos ligados), pior caso de busca (R21).
// =============================================================================

NoLista *inserirLista(NoLista *cab, double area, unsigned int id)
{
    NoLista *n = new NoLista;
    n->bounding_area = area;
    n->sample_id = id;
    n->prox = cab;
    return n;
}
bool buscarLista(NoLista *cab, double area)
{
    for (; cab; cab = cab->prox)
        if (fabs(cab->bounding_area - area) < 1e-6)
            return true;
    return false;
}
NoLista *removerLista(NoLista *cab, unsigned int id)
{
    NoLista *p = cab, *ant = NULL;
    while (p && p->sample_id != id)
    {
        ant = p;
        p = p->prox;
    }
    if (!p)
        return cab;
    if (!ant)
        cab = p->prox;
    else
        ant->prox = p->prox;
    delete p;
    return cab;
}
void liberarLista(NoLista *cab)
{
    while (cab)
    {
        NoLista *t = cab->prox;
        delete cab;
        cab = t;
    }
}

// =============================================================================
//  3. HASH OTIMIZADA (Open Addressing + Tombstone)  [REQUISITO 7: otimizacao]
// -----------------------------------------------------------------------------
//  Versao OTIMIZADA da hash: sondagem linear em array contiguo + tombstone na
//  remocao. Melhor localidade de cache, sem ponteiros (trade-off: mais memoria).
// =============================================================================

TabelaHashOtm *criarHashOtm()
{
    TabelaHashOtm *h = new TabelaHashOtm;
    h->tabela = new HashSlot[TAM_HASH_OTM]();
    return h;
}
void inserirHashOtm(TabelaHashOtm *h, unsigned int id, double area)
{
    int i = id % TAM_HASH_OTM;
    while (h->tabela[i].ocupado)
    {
        if (h->tabela[i].sample_id == id)
        {
            h->tabela[i].bounding_area = area;
            return;
        }
        i = (i + 1) % TAM_HASH_OTM;
    }
    h->tabela[i] = {id, area, true, false};
}
bool buscarHashOtm(TabelaHashOtm *h, unsigned int id, double *area)
{
    int i = id % TAM_HASH_OTM, ini = i;
    while (h->tabela[i].ocupado || h->tabela[i].removido)
    {
        if (h->tabela[i].ocupado && h->tabela[i].sample_id == id)
        {
            *area = h->tabela[i].bounding_area;
            return true;
        }
        i = (i + 1) % TAM_HASH_OTM;
        if (i == ini)
            break;
    }
    return false;
}
void removerHashOtm(TabelaHashOtm *h, unsigned int id)
{
    int i = id % TAM_HASH_OTM, ini = i;
    while (h->tabela[i].ocupado || h->tabela[i].removido)
    {
        if (h->tabela[i].ocupado && h->tabela[i].sample_id == id)
        {
            h->tabela[i].ocupado = false;
            h->tabela[i].removido = true;
            return;
        }
        i = (i + 1) % TAM_HASH_OTM;
        if (i == ini)
            break;
    }
}
void liberarHashOtm(TabelaHashOtm *h)
{
    delete[] h->tabela;
    delete h;
}

// =============================================================================
//  4. ÁRVORE AVL
// -----------------------------------------------------------------------------
//  Arvore de busca auto-balanceada (rotacoes): O(log n) garantido no pior caso
//  e dados ordenados (permite operacoes de faixa, como a poda Z-Score/IQR).
// =============================================================================

int alturaAVL(NoAVL *n) { return n ? n->altura : 0; }
int maxInt(int a, int b) { return a > b ? a : b; }

NoAVL *criarNoAVL(double area, unsigned int id, const char *frame)
{
    NoAVL *n = new NoAVL;
    n->bounding_area = area;
    n->sample_id = id;
    strncpy(n->frame_name, frame, 255);
    n->frame_name[255] = '\0';
    n->esq = n->dir = NULL;
    n->altura = 1;
    return n;
}
int fatorBal(NoAVL *n) { return n ? alturaAVL(n->esq) - alturaAVL(n->dir) : 0; }

NoAVL *rotDir(NoAVL *y)
{
    NoAVL *x = y->esq, *T = x->dir;
    x->dir = y;
    y->esq = T;
    y->altura = maxInt(alturaAVL(y->esq), alturaAVL(y->dir)) + 1;
    x->altura = maxInt(alturaAVL(x->esq), alturaAVL(x->dir)) + 1;
    return x;
}
NoAVL *rotEsq(NoAVL *x)
{
    NoAVL *y = x->dir, *T = y->esq;
    x->dir = T;
    y->esq = x;
    x->altura = maxInt(alturaAVL(x->esq), alturaAVL(x->dir)) + 1;
    y->altura = maxInt(alturaAVL(y->esq), alturaAVL(y->dir)) + 1;
    return y;
}
static NoAVL *balAVL(NoAVL *n)
{
    if (!n)
        return NULL;
    n->altura = 1 + maxInt(alturaAVL(n->esq), alturaAVL(n->dir));
    int b = fatorBal(n);
    if (b > 1 && fatorBal(n->esq) >= 0)
        return rotDir(n);
    if (b > 1 && fatorBal(n->esq) < 0)
    {
        n->esq = rotEsq(n->esq);
        return rotDir(n);
    }
    if (b < -1 && fatorBal(n->dir) <= 0)
        return rotEsq(n);
    if (b < -1 && fatorBal(n->dir) > 0)
    {
        n->dir = rotDir(n->dir);
        return rotEsq(n);
    }
    return n;
}
NoAVL *inserirAVL(NoAVL *n, double area, unsigned int id, const char *frame)
{
    // R2 (Restricao de Memoria): rejeita se a AVL atingiu o teto configurado
    if (g_total_avl_nodes >= g_limite_avl)
    {
        g_avl_rejeitados++;
        return n;
    }
    if (!n)
    {
        g_total_avl_nodes++;
        return criarNoAVL(area, id, frame);
    }
    if (area < n->bounding_area)
        n->esq = inserirAVL(n->esq, area, id, frame);
    else
        n->dir = inserirAVL(n->dir, area, id, frame);
    return balAVL(n);
}
NoAVL *minAVL(NoAVL *n)
{
    while (n->esq)
        n = n->esq;
    return n;
}
NoAVL *removerAVL(NoAVL *r, double area, unsigned int id)
{
    if (!r)
        return NULL;
    if (area < r->bounding_area - 1e-6)
        r->esq = removerAVL(r->esq, area, id);
    else if (area > r->bounding_area + 1e-6)
        r->dir = removerAVL(r->dir, area, id);
    else
    {
        if (r->sample_id != id)
        {
            r->dir = removerAVL(r->dir, area, id);
            return balAVL(r);
        }
        if (!r->esq)
        {
            NoAVL *t = r->dir;
            delete r;
            g_total_avl_nodes--;
            return t;
        }
        if (!r->dir)
        {
            NoAVL *t = r->esq;
            delete r;
            g_total_avl_nodes--;
            return t;
        }
        NoAVL *t = minAVL(r->dir);
        r->bounding_area = t->bounding_area;
        r->sample_id = t->sample_id;
        memcpy(r->frame_name, t->frame_name, sizeof(r->frame_name));
        r->dir = removerAVL(r->dir, t->bounding_area, t->sample_id);
    }
    return balAVL(r);
}
NoAVL *buscarAVL(NoAVL *r, double area)
{
    if (!r)
        return NULL;
    if (fabs(r->bounding_area - area) < 1e-6)
        return r;
    return area < r->bounding_area ? buscarAVL(r->esq, area) : buscarAVL(r->dir, area);
}
void liberarAVL(NoAVL *r)
{
    if (r)
    {
        liberarAVL(r->esq);
        liberarAVL(r->dir);
        delete r;
    }
}

// =============================================================================
//  5. HASH ENCADEAMENTO (Separate Chaining)
// -----------------------------------------------------------------------------
//  Tabela hash classica (encadeamento separado): O(1) medio; base de comparacao
//  para a versao otimizada. Conta colisoes para a metrica do PDF.
// =============================================================================

int funcaoHash(unsigned int id) { return id % TAMANHO_TABELA; }
TabelaHash *criarTabelaHash()
{
    TabelaHash *t = new TabelaHash;
    for (int i = 0; i < TAMANHO_TABELA; i++)
        t->baldes[i] = NULL;
    t->total_colisoes = 0;
    return t;
}
void inserirHash(TabelaHash *t, unsigned int id, double area)
{
    int i = funcaoHash(id);
    if (t->baldes[i])
        t->total_colisoes++;
    for (NoHash *p = t->baldes[i]; p; p = p->prox)
        if (p->sample_id == id)
        {
            p->bounding_area = area;
            return;
        }
    NoHash *n = new NoHash;
    n->sample_id = id;
    n->bounding_area = area;
    n->prox = t->baldes[i];
    t->baldes[i] = n;
}
bool buscarHash(TabelaHash *t, unsigned int id, double *area)
{
    for (NoHash *p = t->baldes[funcaoHash(id)]; p; p = p->prox)
        if (p->sample_id == id)
        {
            *area = p->bounding_area;
            return true;
        }
    return false;
}
bool buscarHashRestrita(TabelaHash *t, unsigned int id, double *area)
{
    if (!g_janela_aberta)
    {
        printf(RED "  [BLOQUEIO R9] CPU ocupada. Busca negada!\n" RST);
        return false;
    }
    return buscarHash(t, id, area);
}
void removerHash(TabelaHash *t, unsigned int id)
{
    int i = funcaoHash(id);
    NoHash *p = t->baldes[i], *ant = NULL;
    while (p && p->sample_id != id)
    {
        ant = p;
        p = p->prox;
    }
    if (!p)
        return;
    if (!ant)
        t->baldes[i] = p->prox;
    else
        ant->prox = p->prox;
    delete p;
}
void liberarTabelaHash(TabelaHash *t)
{
    for (int i = 0; i < TAMANHO_TABELA; i++)
    {
        NoHash *p = t->baldes[i];
        while (p)
        {
            NoHash *nx = p->prox;
            delete p;
            p = nx;
        }
    }
    delete t;
}

// =============================================================================
//  6. KD-TREE (estrutura FORA da ementa #1)
// -----------------------------------------------------------------------------
//  Indexa pontos 2D (cx,cy) alternando a dimensao por nivel. Permite consultas
//  espaciais por frame (range query, vizinho mais proximo, DBSCAN).
// =============================================================================

NoKD *criarNoKD(double pt[], unsigned int id, int frame)
{
    NoKD *n = new NoKD;
    for (int i = 0; i < K; i++)
        n->ponto[i] = pt[i];
    n->sample_id = id;
    n->frame_id = frame;
    n->esq = n->dir = NULL;
    return n;
}
NoKD *inserirKDRec(NoKD *r, double pt[], unsigned int id, int frame, int d)
{
    if (!r)
        return criarNoKD(pt, id, frame);
    int cd = d % K;
    if (pt[cd] < r->ponto[cd])
        r->esq = inserirKDRec(r->esq, pt, id, frame, d + 1);
    else
        r->dir = inserirKDRec(r->dir, pt, id, frame, d + 1);
    return r;
}
NoKD *inserirKD(NoKD *r, double pt[], unsigned int id, int frame) { return inserirKDRec(r, pt, id, frame, 0); }
bool buscarKDRec(NoKD *r, double pt[], int d)
{
    if (!r)
        return false;
    if (fabs(r->ponto[0] - pt[0]) < 1e-6 && fabs(r->ponto[1] - pt[1]) < 1e-6)
        return true;
    int cd = d % K;
    return pt[cd] < r->ponto[cd] ? buscarKDRec(r->esq, pt, d + 1) : buscarKDRec(r->dir, pt, d + 1);
}
bool buscarKD(NoKD *r, double pt[]) { return buscarKDRec(r, pt, 0); }
void liberarKD(NoKD *r)
{
    if (r)
    {
        liberarKD(r->esq);
        liberarKD(r->dir);
        delete r;
    }
}

// =============================================================================
//  7. BLOOM FILTER (estrutura FORA da ementa #2)
// -----------------------------------------------------------------------------
//  Teste de pertinencia probabilistico O(k) com memoria minima (vetor de bits).
//  Sem falso negativo; pode ter falso positivo; nao suporta remocao.
// =============================================================================

unsigned int bf_h1(unsigned int id) { return id * 2654435761u; }
unsigned int bf_h2(unsigned int id)
{
    id ^= id >> 16;
    id *= 0x85ebca6b;
    id ^= id >> 13;
    id *= 0xc2b2ae35;
    id ^= id >> 16;
    return id;
}
unsigned int bf_slot(unsigned int id, int i, int tam)
{
    return (bf_h1(id) + (unsigned int)(i * bf_h2(id))) % (unsigned int)tam;
}

BloomFilter *criarBloom()
{
    BloomFilter *f = new BloomFilter;
    f->tamanho = TAM_BLOOM;
    f->bits = new uint8_t[(TAM_BLOOM + 7) / 8]();
    return f;
}
void inserirBloom(BloomFilter *f, unsigned int id)
{
    for (int i = 0; i < NUM_HASHES; i++)
        SET_BIT(f->bits, bf_slot(id, i, f->tamanho));
}
bool buscarBloom(BloomFilter *f, unsigned int id)
{
    for (int i = 0; i < NUM_HASHES; i++)
        if (!TEST_BIT(f->bits, bf_slot(id, i, f->tamanho)))
            return false;
    return true;
}
void liberarBloom(BloomFilter *f)
{
    delete[] f->bits;
    delete f;
}

// --- Poda da AVL: remove todos os nos com area fora da faixa [mn, mx] ---
NoAVL *podarAVL(NoAVL *r, double mn, double mx)
{
    if (!r)
        return NULL;
    r->esq = podarAVL(r->esq, mn, mx);
    r->dir = podarAVL(r->dir, mn, mx);
    if (r->bounding_area < mn || r->bounding_area > mx)
    {
        if (r->esq && r->dir)
        {
            NoAVL *t = minAVL(r->dir);
            r->bounding_area = t->bounding_area;
            r->sample_id = t->sample_id;
            memcpy(r->frame_name, t->frame_name, sizeof(r->frame_name));
            r->dir = removerAVL(r->dir, t->bounding_area, t->sample_id);
        }
        else
        {
            NoAVL *s = r->esq ? r->esq : r->dir;
            delete r;
            g_total_avl_nodes--;
            return s;
        }
    }
    return balAVL(r);
}

// --- Coleta valores da AVL em in-order para histograma ---
void coletarAVL(NoAVL *r, vector<double> &vals)
{
    if (!r)
        return;
    coletarAVL(r->esq, vals);
    vals.push_back(r->bounding_area);
    coletarAVL(r->dir, vals);
}


// ----- (origem: analise.cpp) -----

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <random>
#include <thread>
#include <chrono>
#include <unordered_map>
using namespace std;

// =============================================================================
//  8. OPERACOES ADICIONAIS (3 exigidas pelo PDF)
//  (a) Estatistica (media/desvio/outliers Z-Score vs IQR); (b) filtragem (poda
//  da AVL); (c) classificacao espacial (densidade por frame na KD-Tree).
// =============================================================================

double calcularMedia(NoLista *cab, int *total)
{
    double soma = 0.0;
    int n = 0;
    for (NoLista *p = cab; p; p = p->prox)
    {
        soma += p->bounding_area;
        n++;
    }
    *total = n;
    return n ? soma / n : 0.0;
}
double calcularDP(NoLista *cab, double media, int total)
{
    if (!cab || total <= 1)
        return 0.0;
    double sq = 0.0;
    for (NoLista *p = cab; p; p = p->prox)
    {
        double d = p->bounding_area - media;
        sq += d * d;
    }
    return sqrt(sq / (total - 1));
}

// --- Estatistica robusta: percentil de um vetor JA ORDENADO (interpolacao linear) ---
static double percentil(const vector<double> &ord, double p)
{
    if (ord.empty())
        return 0.0;
    if (ord.size() == 1)
        return ord[0];
    double idx = p / 100.0 * (ord.size() - 1);
    int lo = (int)idx;
    double frac = idx - lo;
    if (lo + 1 >= (int)ord.size())
        return ord.back();
    return ord[lo] * (1 - frac) + ord[lo + 1] * frac;
}
ResumoEstatistico calcularResumo(vector<double> vals)
{
    ResumoEstatistico r = {};
    if (vals.empty())
        return r;
    sort(vals.begin(), vals.end());

    r.minimo = vals.front();
    r.maximo = vals.back();
    r.q1 = percentil(vals, 25);
    r.mediana = percentil(vals, 50);
    r.q3 = percentil(vals, 75);
    r.iqr = r.q3 - r.q1;

    double soma = 0.0;
    for (double v : vals)
        soma += v;
    r.media = soma / vals.size();

    double sq = 0.0;
    for (double v : vals)
    {
        double d = v - r.media;
        sq += d * d;
    }
    r.desvio = vals.size() > 1 ? sqrt(sq / (vals.size() - 1)) : 0.0;

    // Tukey fences (robusto): Q1 - 1.5*IQR  ate  Q3 + 1.5*IQR
    r.iqr_inf = max(0.0, r.q1 - 1.5 * r.iqr);
    r.iqr_sup = r.q3 + 1.5 * r.iqr;

    // Z-score (sensivel a contaminacao): media +/- 3*sigma
    r.zsc_inf = max(0.0, r.media - 3 * r.desvio);
    r.zsc_sup = r.media + 3 * r.desvio;
    return r;
}

static void painelHistograma(const vector<double> &vals, const char *titulo,
                             const char *cor, bool log_escala)
{
    if (vals.empty())
    {
        printf(YEL "  (sem dados)\n" RST);
        return;
    }
    const int NB = 10;
    double mn = *min_element(vals.begin(), vals.end());
    double mx = *max_element(vals.begin(), vals.end());
    double step = (mx - mn);
    if (step < 1e-9)
    {
        printf(YEL "  Todos os valores sao identicos (%.4f)\n" RST, mn);
        return;
    }
    step /= NB;

    int cnt[10] = {};
    for (double v : vals)
    {
        int b = (int)((v - mn) / step);
        if (b >= NB)
            b = NB - 1;
        cnt[b]++;
    }
    int mx_cnt = *max_element(cnt, cnt + NB);

    printf(CYN BOLD "  %s  (%zu valores)%s\n" RST, titulo, vals.size(),
           log_escala ? "  [escala LOG]" : "  [escala linear]");
    printLinha('-', 57);
    for (int i = 0; i < NB; i++)
    {
        double lo = mn + i * step, hi = lo + step;
        printf(YEL "  [%6.4f-%6.4f] " RST, lo, hi);
        int barras;
        if (log_escala)
            // log1p comprime: balde de 1 ainda aparece, balde de 23k nao estoura
            barras = (cnt[i] > 0) ? (int)(log1p(cnt[i]) / log1p(mx_cnt) * 34) : 0;
        else
            barras = mx_cnt > 0 ? cnt[i] * 34 / mx_cnt : 0;
        printf("%s", cor);
        for (int j = 0; j < barras; j++)
            putchar('#');
        printf(RST " %d\n", cnt[i]);
    }
    printLinha('-', 57);
}

// --- Painel de resumo estatistico (explica a assimetria) ---
static void painelResumo(const ResumoEstatistico &r)
{
    printf(CYN BOLD "  RESUMO ESTATISTICO\n" RST);
    printLinha('-', 57);
    printf("  Minimo : " YEL "%8.4f" RST "    Q1      : " YEL "%8.4f\n" RST, r.minimo, r.q1);
    printf("  Mediana: " YEL "%8.4f" RST "    Media   : " YEL "%8.4f\n" RST, r.mediana, r.media);
    printf("  Q3     : " YEL "%8.4f" RST "    Maximo  : " YEL "%8.4f\n" RST, r.q3, r.maximo);
    printf("  IQR    : " YEL "%8.4f" RST "    Desvio  : " YEL "%8.4f\n" RST, r.iqr, r.desvio);
    printLinha('-', 57);
    // Diagnostico de assimetria: media vs mediana
    if (r.media > r.mediana * 1.5)
        printf(DIM "  Media (%.4f) >> Mediana (%.4f): distribuicao ASSIMETRICA A DIREITA.\n"
                   "  Causa: ~10%% de anomalias R18 puxam a media para cima, mas nao a mediana.\n" RST,
               r.media, r.mediana);
    printLinha('-', 57);
}

void exibirHistograma(NoLista *cab, NoAVL *avl)
{
    if (!cab)
    {
        printf(RED "  Sem dados. Ingira o dataset primeiro.\n" RST);
        return;
    }

    vector<double> vals_lista;
    for (NoLista *p = cab; p; p = p->prox)
        vals_lista.push_back(p->bounding_area);

    // Resumo estatistico primeiro: explica a assimetria
    printf("\n");
    ResumoEstatistico r = calcularResumo(vals_lista);
    painelResumo(r);

    // Painel 1: linear (mostra o pico dominante = dados reais)
    printf("\n");
    painelHistograma(vals_lista, "DADOS BRUTOS - linear", GRN, false);

    // Painel 2: LOG (revela a estrutura escondida = cauda de anomalias R18)
    printf("\n");
    painelHistograma(vals_lista, "DADOS BRUTOS - logaritmico", CYN, true);

    // Painel 3: AVL atual (reflete poda se ja executada)
    vector<double> vals_avl;
    coletarAVL(avl, vals_avl);
    printf("\n");
    if (vals_avl.size() == vals_lista.size())
        painelHistograma(vals_avl, "AVL (poda nao executada ainda)", YEL, true);
    else
        painelHistograma(vals_avl, "AVL (pos-poda - outliers removidos)", BLU, true);

    printf(DIM "  O painel LOG revela que os baldes 'vazios' tem dezenas a centenas\n"
               "  de anomalias R18. Rode a opcao 3 para filtra-las e compare a AVL.\n" RST);
}
// =============================================================================
//  ANALISE ESPACIAL VIA KD-TREE
//  cx,cy sao a posicao da deteccao no quadro da imagem (por frame), nao a
//  geografia do campo. (A) distribuicao, (B) vizinho, (C) produtividade/cluster.
// =============================================================================

// --- Range query: conta pontos dentro de um retangulo, podando subarvores ---
void contarZona(NoKD *r, double x0, double x1, double y0, double y1, int d, int *c)
{
    if (!r)
        return;
    double x = r->ponto[0], y = r->ponto[1];
    if (x >= x0 && x < x1 && y >= y0 && y < y1)
        (*c)++;
    int cd = d % K;
    if (cd == 0)
    {
        if (x0 < x)
            contarZona(r->esq, x0, x1, y0, y1, d + 1, c);
        if (x1 >= x)
            contarZona(r->dir, x0, x1, y0, y1, d + 1, c);
    }
    else
    {
        if (y0 < y)
            contarZona(r->esq, x0, x1, y0, y1, d + 1, c);
        if (y1 >= y)
            contarZona(r->dir, x0, x1, y0, y1, d + 1, c);
    }
}

// (A) Distribuicao espacial das deteccoes NO QUADRO (nao no campo)
void analisarDistribuicaoEspacial(NoKD *kd)
{
    int q[4] = {};
    contarZona(kd, 0.0, 0.5, 0.5, 1.0, 0, &q[0]);
    contarZona(kd, 0.5, 1.0, 0.5, 1.0, 0, &q[1]);
    contarZona(kd, 0.0, 0.5, 0.0, 0.5, 0, &q[2]);
    contarZona(kd, 0.5, 1.0, 0.0, 0.5, 0, &q[3]);

    int total = q[0] + q[1] + q[2] + q[3];
    if (!total)
    {
        printf(YEL "  Nenhum ponto em [0,1]. Verifique coordenadas.\n" RST);
        return;
    }
    double med = total / 4.0;
    const char *nomes[4] = {"Quadro Sup-Esquerda", "Quadro Sup-Direita",
                            "Quadro Inf-Esquerda", "Quadro Inf-Direita"};

    printf("\n" CYN BOLD "  (A) DISTRIBUICAO DAS DETECCOES NO QUADRO (KD-Tree range query)\n" RST);
    printLinha('-', 57);
    int mx = *max_element(q, q + 4);
    for (int i = 0; i < 4; i++)
    {
        const char *cor = (q[i] >= med * 1.2) ? GRN : (q[i] <= med * 0.8) ? RED
                                                                          : YEL;
        printf("  %-20s ", nomes[i]);
        printBarra(q[i], mx, 18, cor);
        printf(" %5d", q[i]);
        if (q[i] >= med * 1.2)
            printf(GRN " [concentracao ALTA]\n" RST);
        else if (q[i] <= med * 0.8)
            printf(RED " [concentracao BAIXA]\n" RST);
        else
            printf(YEL " [normal]\n" RST);
    }
    printLinha('-', 57);
    printf(DIM "  Mede a regiao DO QUADRO onde as uvas aparecem (vies de camera),\n"
               "  nao a geografia do campo - cx,cy sao coordenadas da imagem.\n" RST);
}

// --- (B) Busca do VIZINHO MAIS PROXIMO (nearest-neighbor) na KD-Tree ---------
static double distSqKD(const double a[], const double b[])
{
    double dx = a[0] - b[0], dy = a[1] - b[1];
    return dx * dx + dy * dy;
}
void nnKD(NoKD *r, const double alvo[], int depth,
          double *melhorDist, unsigned int *melhorId, unsigned int exclId)
{
    if (!r)
        return;
    if (r->sample_id != exclId)
    {
        double dd = distSqKD(r->ponto, alvo);
        if (dd < *melhorDist)
        {
            *melhorDist = dd;
            *melhorId = r->sample_id;
        }
    }
    int cd = depth % K;
    double diff = alvo[cd] - r->ponto[cd];
    NoKD *perto = diff < 0 ? r->esq : r->dir;
    NoKD *longe = diff < 0 ? r->dir : r->esq;
    nnKD(perto, alvo, depth + 1, melhorDist, melhorId, exclId);
    // so explora o outro lado se ainda houver chance de algo mais perto la
    if (diff * diff < *melhorDist)
        nnKD(longe, alvo, depth + 1, melhorDist, melhorId, exclId);
}
static void coletarPontosKD(NoKD *r, vector<PontoKD> &v)
{
    if (!r)
        return;
    coletarPontosKD(r->esq, v);
    v.push_back({r->ponto[0], r->ponto[1], r->sample_id, r->frame_id});
    coletarPontosKD(r->dir, v);
}

// (B) Vizinho mais proximo: dado um ponto, acha a deteccao mais proxima O(log n).
void analisarVizinhanca(NoKD *kd)
{
    printf("\n" CYN BOLD "  (B) BUSCA DO VIZINHO MAIS PROXIMO (KD-Tree NN)\n" RST);
    printLinha('-', 57);
    printf("  Dado um ponto (x,y) no quadro, acha a deteccao mais proxima:\n");

    struct Q
    {
        double x, y;
        const char *nome;
    };
    Q consultas[] = {
        {0.50, 0.50, "Centro"},
        {0.10, 0.10, "Sup-esquerda"},
        {0.90, 0.90, "Inf-direita"},
        {0.90, 0.10, "Sup-direita"},
    };
    for (const Q &q : consultas)
    {
        double pt[2] = {q.x, q.y};
        double bd = 1e18;
        unsigned int bid = 0;
        nnKD(kd, pt, 0, &bd, &bid, 0);
        if (bid)
            printf("    (%.2f, %.2f) %-13s -> " YEL "ID %u" RST " a dist %.4f\n",
                   q.x, q.y, q.nome, bid, sqrt(bd));
    }
    printLinha('-', 57);
    printf(DIM "  A poda da KD-Tree torna a busca O(log n) medio, contra O(n) de\n"
               "  comparar com todas as deteccoes uma a uma.\n" RST);
}

// Range query por raio: coleta ids dentro de um circulo (filtrando por frame).
static void rangeKD(NoKD *r, const double alvo[], double eps, int frame, int depth,
                    vector<unsigned int> &out)
{
    if (!r)
        return;
    if ((frame < 0 || r->frame_id == frame) && distSqKD(r->ponto, alvo) <= eps * eps)
        out.push_back(r->sample_id);
    int cd = depth % K;
    double diff = alvo[cd] - r->ponto[cd];
    NoKD *perto = diff < 0 ? r->esq : r->dir;
    NoKD *longe = diff < 0 ? r->dir : r->esq;
    rangeKD(perto, alvo, eps, frame, depth + 1, out);
    if (diff * diff <= eps * eps)
        rangeKD(longe, alvo, eps, frame, depth + 1, out);
}

// (C) Produtividade por frame: conta cachos por frame e classifica na faixa
// agronomica (poucos = subproducao, muitos = sobrecarga).
void produtividadePorFrame(const vector<PontoKD> &pts)
{
    const int IDEAL_MIN = 9;
    const int IDEAL_MAX = 20;

    unordered_map<int, int> cachosPorFrame;
    for (const PontoKD &p : pts)
        cachosPorFrame[p.frame]++;

    int sub = 0, ideal = 0, sobre = 0, maxc = 0;
    double soma = 0;
    for (auto &kv : cachosPorFrame)
    {
        int c = kv.second;
        if (c < IDEAL_MIN)
            sub++;
        else if (c <= IDEAL_MAX)
            ideal++;
        else
            sobre++;
        soma += c;
        maxc = max(maxc, c);
    }
    int nframes = (int)cachosPorFrame.size();
    if (!nframes)
        return;

    printf("\n" CYN BOLD "  (C) PRODUTIVIDADE POR FRAME (cachos por imagem/planta)\n" RST);
    printLinha('-', 57);
    printf("  Frames analisados:  " YEL "%d" RST "   | media: " YEL "%.1f" RST
           " cachos/frame  | max: " YEL "%d\n" RST,
           nframes, soma / nframes, maxc);
    printLinha('-', 57);
    printf("  Classificacao agronomica (faixa ideal %d-%d cachos/frame):\n", IDEAL_MIN, IDEAL_MAX);
    int mx = max({sub, ideal, sobre, 1});
    printf("    " RED "Subproducao (< %2d): " RST, IDEAL_MIN);
    printBarra(sub, mx, 22, RED);
    printf(" %d frames\n", sub);
    printf("    " GRN "Ideal (%2d-%2d):      " RST, IDEAL_MIN, IDEAL_MAX);
    printBarra(ideal, mx, 22, GRN);
    printf(" %d frames\n", ideal);
    printf("    " YEL "Sobrecarga (> %2d):  " RST, IDEAL_MAX);
    printBarra(sobre, mx, 22, YEL);
    printf(" %d frames -> raleio\n", sobre);
    printLinha('-', 57);
    printf(DIM "  Agora cada frame tem coordenadas proprias (frame_id recuperado dos\n"
               "  nomes dos arquivos do Roboflow), entao cachos/frame e uma metrica\n"
               "  agronomica VALIDA - nao mistura imagens diferentes.\n" RST);
}

// (D) DBSCAN por frame: agrupa cachos adjacentes do mesmo frame (KD-Tree faz as
// vizinhancas em O(log n) -> ~O(n log n), contra O(n^2) da forca bruta).
void clusterizarPorFrame(NoKD *kd, const vector<PontoKD> &pts,
                         const unordered_map<unsigned int, int> &idx)
{
    const int MINPTS = 3;
    const double EPS = 0.10;

    int m = (int)pts.size();
    vector<int> label(m, -1);
    int cid = 0;
    for (int i = 0; i < m; i++)
    {
        if (label[i] != -1)
            continue;
        double alvo[2] = {pts[i].x, pts[i].y};
        vector<unsigned int> viz;
        rangeKD(kd, alvo, EPS, pts[i].frame, 0, viz);
        if ((int)viz.size() < MINPTS)
        {
            label[i] = 0;
            continue;
        }
        cid++;
        label[i] = cid;
        vector<int> fila;
        for (unsigned int id : viz)
        {
            int j = idx.at(id);
            if (j != i)
                fila.push_back(j);
        }
        for (size_t s = 0; s < fila.size(); s++)
        {
            int j = fila[s];
            if (label[j] == 0)
                label[j] = cid;
            if (label[j] != -1)
                continue;
            label[j] = cid;
            double a2[2] = {pts[j].x, pts[j].y};
            vector<unsigned int> viz2;
            rangeKD(kd, a2, EPS, pts[j].frame, 0, viz2);
            if ((int)viz2.size() >= MINPTS)
                for (unsigned int id2 : viz2)
                    fila.push_back(idx.at(id2));
        }
    }
    int isolados = 0;
    for (int i = 0; i < m; i++)
        if (label[i] <= 0)
            isolados++;

    printf("\n" CYN BOLD "  (D) GRUPOS DE CACHOS COLADOS POR FRAME (DBSCAN c/ KD-Tree)\n" RST);
    printLinha('-', 57);
    printf("  eps=%.2f minPts=%d (so conecta cachos do MESMO frame)\n", EPS, MINPTS);
    printf("  Grupos de cachos adjacentes: " YEL "%d" RST "   | cachos isolados: " YEL "%d\n" RST,
           cid, isolados);
    printLinha('-', 57);
    printf(DIM "  Filtrando por frame, o DBSCAN agrupa so cachos realmente vizinhos\n"
               "  na MESMA imagem (cachos colados sao reais, nao erro). A KD-Tree faz\n"
               "  cada consulta de vizinhanca em O(log n) -> ~O(n log n) total.\n" RST);
}

// Relatorio espacial completo (A + B + C + D)
void analiseEspacialKD(NoKD *kd)
{
    analisarDistribuicaoEspacial(kd);
    analisarVizinhanca(kd);

    vector<PontoKD> pts;
    coletarPontosKD(kd, pts);
    if (pts.empty())
        return;
    unordered_map<unsigned int, int> idx;
    idx.reserve(pts.size() * 2);
    for (int i = 0; i < (int)pts.size(); i++)
        idx[pts[i].id] = i;

    produtividadePorFrame(pts);
    clusterizarPorFrame(kd, pts, idx);
}

// =============================================================================
//  9. GERADOR DE DATASET SINTÉTICO
// =============================================================================

void gerarDataset(const string &nome, int n)
{
    printf("  Gerando %d amostras sinteticas em '%s'...\n", n, nome.c_str());
    ofstream f(nome);
    if (!f.is_open())
    {
        printf(RED "  [ERRO] Falha ao criar arquivo.\n" RST);
        return;
    }
    // Schema novo: sample_id, frame_id, cx, cy, width, height, area
    f << "sample_id,frame_id,cx,cy,width,height,area\n";
    mt19937 gen(2026);
    uniform_real_distribution<double> dxy(0.05, 0.95), dwh(0.01, 0.15);
    uniform_int_distribution<int> dCachos(5, 20);
    normal_distribution<double> ruido(0.0, 0.003);
    int sid = 0, frame = 0;
    while (sid < n)
    {
        int cachos = dCachos(gen);
        for (int c = 0; c < cachos && sid < n; c++)
        {
            double cx = dxy(gen), cy = dxy(gen), w = dwh(gen), h = dwh(gen);
            double area = max(0.0001, w * h + ruido(gen));
            f << sid << "," << frame << "," << cx << "," << cy << "," << w << "," << h << "," << area << "\n";
            sid++;
        }
        frame++;
    }
    f.close();
    printf(GRN "  [OK] Dataset gerado com sucesso.\n" RST);
}

// =============================================================================
//  10. INGESTÃO DO CSV
// =============================================================================

void carregarDataset(NoLista **lista, TabelaHashOtm *hotm, NoAVL **avl,
                     TabelaHash *hash, NoKD **kd, BloomFilter *bloom)
{
    string nome = "dataset_uvas_consolidado.csv";
    ifstream arq(nome);
    if (!arq.is_open())
    {
        printf(YEL "\n  [AVISO] '%s' nao encontrado.\n" RST, nome.c_str());
        printf("  Gerar dataset sintetico com 500.000 amostras? (s/n): ");
        char r;
        cin >> r;
        if (r == 's' || r == 'S')
        {
            gerarDataset(nome, 500000);
            arq.open(nome);
        }
        if (!arq.is_open())
            return;
    }

    // conta linhas para barra de progresso
    unsigned int total_linhas = 0;
    {
        string ln;
        getline(arq, ln);
        while (getline(arq, ln))
            total_linhas++;
    }
    arq.clear();
    arq.seekg(0);
    string cabecalho;
    getline(arq, cabecalho);

    printf("\n" CYN "  Ingerindo '%s' (%u registros)...\n" RST, nome.c_str(), total_linhas);

    mt19937 gen(42);
    uniform_real_distribution<double> dNoise(0.40, 0.95);
    uniform_int_distribution<int> dDefeito(1, 100);

    unsigned int cnt = 0;
    string linha;
    while (getline(arq, linha))
    {
        if (!linha.empty() && linha.back() == '\r')
            linha.pop_back();
        if (linha.empty())
            continue;
        stringstream ss(linha);
        string cel;
        vector<string> cols;
        while (getline(ss, cel, ','))
            cols.push_back(cel);
        // Schema novo: sample_id, frame_id, cx, cy, width, height, area (7 colunas)
        if (cols.size() < 7)
            continue;
        try
        {
            unsigned int id = stoul(cols[0]);
            int frame = stoi(cols[1]);
            double cx = stod(cols[2]), cy = stod(cols[3]), area = stod(cols[6]);
            if (dDefeito(gen) <= 10)
                area = dNoise(gen);
            if (cnt % 1000 == 0)
                this_thread::sleep_for(chrono::milliseconds(1));

            *lista = inserirLista(*lista, area, id);
            inserirHashOtm(hotm, id, area);
            *avl = inserirAVL(*avl, area, id, "frame_real.png");
            inserirHash(hash, id, area);
            double pt[2] = {cx, cy};
            *kd = inserirKD(*kd, pt, id, frame);
            inserirBloom(bloom, id);
            cnt++;
            if (cnt % 10000 == 0)
                printProgresso(cnt, total_linhas);
        }
        catch (...)
        {
            continue;
        }
    }
    arq.close();
    printProgresso(cnt, total_linhas);
    printf("\n" GRN "\n  [OK] %u amostras carregadas.\n" RST, cnt);
}


// ----- (origem: benchmark.cpp) -----

#include <cstdio>
#include <chrono>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
using namespace std;

// =============================================================================
//  11. BENCHMARK COMPLETO
//  Metodologia: (1) so compara estruturas comparaveis (Grupo A = busca por
//  chave; Bloom no Grupo B); (2) alvo variado entre 64 pontos (evita O(1)
//  enganoso e hoisting); (3) repeticoes adaptativas p/ superar o timer (~0.6ms).
// =============================================================================

// Exibe uma linha de benchmark com barra proporcional -usa nanosegundos internamente
static void barraComparativa(const char *nome, double ns, double mx, int larg = 28)
{
    printf("  %-28s ", nome);
    printBarra(ns, mx, larg);
    if (ns < 1000.0)
        printf(" %8.2f ns/op\n", ns);
    else
        printf(" %8.2f us/op\n", ns / 1000.0);
}

// --- Medicao robusta de remocao em hash ---
static NoHash *g_br_baldes[TAMANHO_TABELA];
static HashSlot g_br_slots[TAM_HASH_OTM];

static double medirRemChain(const vector<unsigned int> &ids, const vector<double> &ar)
{
    int n = (int)ids.size();
    for (int b = 0; b < TAMANHO_TABELA; b++)
        g_br_baldes[b] = NULL;
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int b = id % TAMANHO_TABELA;
        NoHash *x = new NoHash;
        x->sample_id = id;
        x->bounding_area = ar[i];
        x->prox = g_br_baldes[b];
        g_br_baldes[b] = x;
    }
    asm volatile("" ::: "memory");
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int b = id % TAMANHO_TABELA;
        NoHash *p = g_br_baldes[b], *ant = NULL;
        while (p && p->sample_id != id)
        {
            ant = p;
            p = p->prox;
        }
        if (p)
        {
            if (!ant)
                g_br_baldes[b] = p->prox;
            else
                ant->prox = p->prox;
            delete p;
        }
    }
    long long tot = chrono::duration_cast<chrono::nanoseconds>(
                        chrono::high_resolution_clock::now() - t1)
                        .count();
    return (double)tot / n;
}

static double medirRemOpen(const vector<unsigned int> &ids, const vector<double> &ar)
{
    int n = (int)ids.size();
    for (int s = 0; s < TAM_HASH_OTM; s++)
    {
        g_br_slots[s].ocupado = false;
        g_br_slots[s].removido = false;
    }
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int s = id % TAM_HASH_OTM;
        while (g_br_slots[s].ocupado)
            s = (s + 1) % TAM_HASH_OTM;
        g_br_slots[s] = {id, ar[i], true, false};
    }
    asm volatile("" ::: "memory");
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++)
    {
        unsigned int id = ids[i];
        int s = id % TAM_HASH_OTM, ini = s;
        while (g_br_slots[s].ocupado || g_br_slots[s].removido)
        {
            if (g_br_slots[s].ocupado && g_br_slots[s].sample_id == id)
            {
                g_br_slots[s].ocupado = false;
                g_br_slots[s].removido = true;
                break;
            }
            s = (s + 1) % TAM_HASH_OTM;
            if (s == ini)
                break;
        }
    }
    long long tot = chrono::duration_cast<chrono::nanoseconds>(
                        chrono::high_resolution_clock::now() - t1)
                        .count();
    // Varre o array inteiro contando tombstones -> observa TODAS as remocoes
    long long chk = 0;
    for (int s = 0; s < TAM_HASH_OTM; s++)
        if (g_br_slots[s].removido)
            chk++;
    g_sink_d += chk;
    return (double)tot / n;
}

void executarBenchmark(NoLista *lista, TabelaHashOtm *hotm, NoAVL *avl,
                       TabelaHash *hash, BloomFilter *bloom, unsigned int id_ref)
{
    // Conta a lista e usa o elemento do meio como ID de referencia no cabecalho.
    int n_lista = 0;
    for (NoLista *p = lista; p; p = p->prox)
        n_lista++;
    NoLista *meio = lista;
    for (int i = 0; i < n_lista / 2 && meio->prox; i++)
        meio = meio->prox;
    id_ref = meio->sample_id;

    // Coleta 64 alvos (id, area) distribuidos pelo dataset. Variar o alvo a cada
    const int NP = 64;
    unsigned int pid[NP];
    double parea[NP];
    {
        NoLista *p = lista;
        int passo = max(1, n_lista / NP);
        for (int k = 0; k < NP; k++)
        {
            // avanca 'passo' nos (ou reinicia se a lista acabar)
            for (int s = 0; s < passo && p && p->prox; s++)
                p = p->prox;
            if (!p)
                p = lista;
            pid[k] = p->sample_id;
            parea[k] = p->bounding_area;
        }
    }

    // Mede ns/op com reps ADAPTATIVO. O alvo cicla via volatile (sem hoisting).
    auto medirR = [&](int reps, auto fn) -> double
    {
        volatile int ctr = 0;
        auto t1 = chrono::high_resolution_clock::now();
        for (int i = 0; i < reps; i++)
        {
            int k = (ctr++) % NP;
            fn(pid[k], parea[k]);
        }
        long long tot = chrono::duration_cast<chrono::nanoseconds>(
                            chrono::high_resolution_clock::now() - t1)
                            .count();
        return (double)tot / reps;
    };
    const int REP_LENTO = 3000, REP_RAPIDO = 2000000;
    // cronoNs(): cronometra um bloco e divide pelo nº de operacoes -> ns/op.
    auto cronoNs = [&](int divisor, auto fn) -> double
    {
        auto t1 = chrono::high_resolution_clock::now();
        fn();
        asm volatile("" ::: "memory");
        return (double)chrono::duration_cast<chrono::nanoseconds>(
                   chrono::high_resolution_clock::now() - t1)
                   .count() /
               divisor;
    };

    double ns_henc = medirR(REP_RAPIDO, [&](unsigned int id, double)
                            { double a=0; g_sink_b = buscarHash(hash, id, &a); g_sink_d=a; });
    double ns_hotm = medirR(REP_RAPIDO, [&](unsigned int id, double)
                            { double a=0; g_sink_b = buscarHashOtm(hotm, id, &a); g_sink_d=a; });
    double ns_avl = medirR(REP_RAPIDO, [&](unsigned int, double ar)
                           { g_sink_b = (buscarAVL(avl, ar) != NULL); });
    double ns_lst = medirR(REP_LENTO, [&](unsigned int, double ar)
                           { g_sink_b = buscarLista(lista, ar); });
    double ns_blm = medirR(REP_RAPIDO, [&](unsigned int id, double)
                           { g_sink_b = buscarBloom(bloom, id); });

    printf("\n" CYN BOLD "  BENCHMARK COMPLETO - todas as estruturas implementadas\n" RST);
    printf(DIM "  (cada estrutura medida isoladamente; valores em ns por operacao)\n" RST);
    printLinha('=', 57);

    // GRUPO A: estruturas que respondem a MESMA pergunta ("achar chave X").
    double mx_a = max({ns_henc, ns_hotm, ns_avl, ns_lst, 0.01});
    printf(BOLD "\n  [ GRUPO A - BUSCA POR CHAVE (comparavel) ]" RST);
    printf(DIM " alvo: ID %u\n" RST, id_ref);
    printLinha('-', 57);
    barraComparativa("Hash Encadeamento O(1):", ns_henc, mx_a);
    barraComparativa("Hash Open Addr. (otim.):", ns_hotm, mx_a);
    barraComparativa("AVL O(log n):", ns_avl, mx_a);
    barraComparativa("Lista O(n) [R21 degradado]:", ns_lst, mx_a);
    printf(DIM "  Todas resolvem 'encontrar o elemento de chave X' -> comparacao valida.\n" RST);

    // GRUPO B: estrutura de PROPOSITO DISTINTO. NAO comparar na mesma regua.
    printf(BOLD "\n  [ GRUPO B - TESTE DE PERTINENCIA (NAO comparavel) ]\n" RST);
    printLinha('-', 57);
    barraComparativa("Bloom Filter O(k):", ns_blm, max(ns_blm, 0.01));
    printf(DIM "  Bloom responde 'X PODE estar no conjunto?' (com falsos positivos),\n"
               "  nao 'onde esta X?'. Por isso NAO entra no ranking do Grupo A.\n"
               "  Idem KD-Tree (busca espacial 2D) - ver opcao 2.\n" RST);

    // ===== INSERCAO: TODAS as estruturas (em copias TEMPORARIAS isoladas) =====
    const int N_INS = 300000;
    vector<double> areas(N_INS);
    {
        mt19937 g(99);
        uniform_real_distribution<double> d(0.001, 0.999);
        for (double &a : areas)
            a = d(g);
    }

    NoLista *t_lst = NULL;
    double ins_lst = cronoNs(N_INS, [&]
                             { for (int i=0;i<N_INS;i++) t_lst = inserirLista(t_lst, areas[i], (unsigned)i); });
    liberarLista(t_lst);

    double ins_avl = 0;
    {
        EstadoAVLGuard g;
        g_total_avl_nodes = 0;
        g_limite_avl = N_INS + 10;
        NoAVL *t = NULL;
        ins_avl = cronoNs(N_INS, [&]
                          { for (int i=0;i<N_INS;i++) t = inserirAVL(t, areas[i], (unsigned)i, "bench"); });
        liberarAVL(t);
    }

    TabelaHash *t_hc = criarTabelaHash();
    double ins_hc = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirHash(t_hc, (unsigned)i, areas[i]); });
    liberarTabelaHash(t_hc);

    TabelaHashOtm *t_ho = criarHashOtm();
    double ins_ho = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirHashOtm(t_ho, (unsigned)i, areas[i]); });
    liberarHashOtm(t_ho);

    NoKD *t_kd = NULL;
    double ins_kd = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++){ double pt[2]={areas[i], areas[(i*7)%N_INS]}; t_kd = inserirKD(t_kd, pt, (unsigned)i, 0); } });
    liberarKD(t_kd);

    BloomFilter *t_bl = criarBloom();
    double ins_bl = cronoNs(N_INS, [&]
                            { for (int i=0;i<N_INS;i++) inserirBloom(t_bl, (unsigned)i); });
    liberarBloom(t_bl);

    double mx_i = max({ins_lst, ins_avl, ins_hc, ins_ho, ins_kd, ins_bl, 0.01});
    printf(BOLD "\n  [ INSERCAO - todas as estruturas (%d elementos cada) ]\n" RST, N_INS);
    printLinha('-', 57);
    barraComparativa("Lista (prepend O(1)):", ins_lst, mx_i);
    barraComparativa("AVL O(log n):", ins_avl, mx_i);
    barraComparativa("Hash Encadeamento O(1):", ins_hc, mx_i);
    barraComparativa("Hash Open Addr. O(1):", ins_ho, mx_i);
    barraComparativa("KD-Tree O(log n):", ins_kd, mx_i);
    barraComparativa("Bloom Filter O(k):", ins_bl, mx_i);

    // ===== REMOCAO: estruturas que suportam (Bloom e KD = N/A justificado) =====
    const int N_REM = 10000;
    double rem_lst = 0, rem_avl = 0, rem_hc = 0, rem_ho = 0;
    // IDs pseudo-aleatorios: o compilador nao consegue simular as estruturas em
    vector<unsigned int> ids(N_REM);
    {
        mt19937 g(123);
        for (auto &x : ids)
            x = g();
    }

    {
        NoLista *t = NULL;
        for (int i = 0; i < N_REM; i++)
            t = inserirLista(t, areas[i], ids[i]);
        rem_lst = cronoNs(N_REM, [&]
                          { for (int i=0;i<N_REM;i++) t = removerLista(t, ids[i]); });
        liberarLista(t);
    }

    {
        EstadoAVLGuard g;
        g_total_avl_nodes = 0;
        g_limite_avl = N_REM + 10;
        NoAVL *t = NULL;
        for (int i = 0; i < N_REM; i++)
            t = inserirAVL(t, areas[i], ids[i], "bench");
        rem_avl = cronoNs(N_REM, [&]
                          { for (int i=0;i<N_REM;i++) t = removerAVL(t, areas[i], ids[i]); });
        liberarAVL(t);
    }

    // Hash: medida em buffers estaticos (medirRemChain/medirRemOpen). Usa MUITOS
    const int N_HREM = 300000;
    vector<unsigned int> hids(N_HREM);
    vector<double> hareas(N_HREM);
    {
        mt19937 g(123);
        uniform_real_distribution<double> d(0.001, 0.999);
        for (int i = 0; i < N_HREM; i++)
        {
            hids[i] = g();
            hareas[i] = d(g);
        }
    }
    rem_hc = medirRemChain(hids, hareas);
    rem_ho = medirRemOpen(hids, hareas);

    double mx_r = max({rem_lst, rem_avl, rem_hc, rem_ho, 0.01});
    printf(BOLD "\n  [ REMOCAO - estruturas que suportam ]\n" RST);
    printLinha('-', 57);
    barraComparativa("Lista O(n):", rem_lst, mx_r);
    barraComparativa("AVL O(log n):", rem_avl, mx_r);
    barraComparativa("Hash Encadeamento O(1):", rem_hc, mx_r);
    barraComparativa("Hash Open Addr. O(1):", rem_ho, mx_r);
    printf("  %-28s " YEL "N/A" RST DIM " - Bloom nao remove (apagar bits afetaria outros ids)\n" RST, "Bloom Filter:");
    printf("  %-28s " YEL "N/A" RST DIM " - KD-Tree: remocao exige reconstruir a subarvore\n" RST, "KD-Tree:");

    // MEMÓRIA ESTIMADA
    long long mem_henc = (long long)TAMANHO_TABELA * sizeof(NoHash *) + sizeof(TabelaHash);
    long long mem_hotm = (long long)TAM_HASH_OTM * sizeof(HashSlot);
    long long mem_blm = (TAM_BLOOM + 7) / 8;
    long long mem_avl = (long long)g_total_avl_nodes * sizeof(NoAVL);
    long long mx_mem = max({mem_henc, mem_hotm, mem_blm, mem_avl});

    printf(BOLD "\n  [ USO DE MEMORIA POR ESTRUTURA (estimado via sizeof) ]\n" RST);
    printLinha('-', 57);
    auto barraKB = [&](const char *n, long long bytes)
    {
        printf("  %-28s ", n);
        printBarra((double)bytes, (double)mx_mem, 28, BLU);
        printf(" %6lld KB\n", bytes / 1024);
    };
    barraKB("Hash Encadeamento (tabela):", mem_henc);
    barraKB("Hash Open Addressing:", mem_hotm);
    barraKB("Bloom Filter:", mem_blm);
    barraKB("AVL (nos alocados):", mem_avl);

    // MEMORIA REAL do processo (medida pelo SO, nao estimada)
    long long rss = memoriaRSS_KB();
    printf(BOLD "\n  [ MEMORIA RAM REAL DO PROCESSO (medida pelo SO) ]\n" RST);
    printLinha('-', 57);
    if (rss >= 0)
    {
        printf("  Working Set (RSS) atual: " GRN "%lld KB" RST " (~%.1f MB)\n",
               rss, rss / 1024.0);
        long long soma_est = (mem_henc + mem_hotm + mem_blm + mem_avl) / 1024;
        printf("  Soma estimada (sizeof):  " YEL "%lld KB\n" RST, soma_est);
        printf(DIM "  RSS conta so paginas RESIDENTES. As tabelas hash sao zeradas:\n"
                   "  slots vazios (nunca escritos) nao entram no working set por\n"
                   "  causa da alocacao preguicosa do SO -> RSS pode ser < estimativa.\n" RST);
    }
    else
    {
        printf(YEL "  (medicao de RSS indisponivel nesta plataforma)\n" RST);
    }

    // COLISÕES + fator de carga explicado
    long long n_elem = g_total_avl_nodes;
    double fator_carga = (TAMANHO_TABELA > 0) ? (double)n_elem / TAMANHO_TABELA * 100.0 : 0.0;

    // Demo de colisao com tabela pequena (mostra o fenomeno em escala)
    const int DEMO_TAM = 97;
    long long colisoes_demo = 0;
    {
        vector<NoHash *> demo(DEMO_TAM, nullptr);
        mt19937 g(7);
        uniform_int_distribution<int> dID(1, 500);
        for (int i = 0; i < 200; i++)
        {
            int id = dID(g);
            int slot = id % DEMO_TAM;
            if (demo[slot])
                colisoes_demo++;
            // nao aloca de verdade, so conta
        }
    }

    printf(BOLD "\n  [ COLISOES - Hash Encadeamento ]\n" RST);
    printLinha('-', 57);
    printf("  Tabela principal (%d slots, %lld elementos):\n", TAMANHO_TABELA, n_elem);
    printf("  Fator de carga:     " YEL "%.2f%%\n" RST, fator_carga);
    printf("  Colisoes reais:     " YEL "%lld\n" RST, hash->total_colisoes);
    if (fator_carga < 10.0)
        printf(DIM "  (Fator de carga baixo -> colisoes raras. Normal para este dataset.)\n" RST);
    printf("\n  Demo colisao (tabela de %d slots, 200 insercoes):\n", DEMO_TAM);
    printf("  Colisoes esperadas: " YEL "%lld (fator de carga ~206%%)\n" RST, colisoes_demo);
    printf(DIM "  Quanto maior o fator de carga, mais colisoes e mais lenta a busca.\n" RST);

    // LATÊNCIA MÉDIA DE OPERAÇÕES COMBINADAS (metrica exigida pelo PDF)
    const int N_CICLOS = 50000;
    auto t_ciclo1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N_CICLOS; i++)
    {
        unsigned int id = 7000000 + i;
        double a = 0.033;
        inserirHash(hash, id, a);
        g_sink_b = buscarHash(hash, id, &a);
        removerHash(hash, id);
    }
    long long lat_henc = chrono::duration_cast<chrono::nanoseconds>(
                             chrono::high_resolution_clock::now() - t_ciclo1)
                             .count() /
                         N_CICLOS;

    auto t_ciclo2 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N_CICLOS; i++)
    {
        unsigned int id = 7000000 + i;
        double a = 0.033;
        inserirHashOtm(hotm, id, a);
        g_sink_b = buscarHashOtm(hotm, id, &a);
        removerHashOtm(hotm, id);
    }
    long long lat_hotm = chrono::duration_cast<chrono::nanoseconds>(
                             chrono::high_resolution_clock::now() - t_ciclo2)
                             .count() /
                         N_CICLOS;

    printf(BOLD "\n  [ LATENCIA MEDIA - operacoes combinadas (ins+busca+rem) ]\n" RST);
    printLinha('-', 57);
    long long mx_lat = max(lat_henc, lat_hotm);
    if (mx_lat == 0)
        mx_lat = 1;
    printf("  %-28s ", "Hash Encadeamento:");
    printBarra(lat_henc, mx_lat, 28, CYN);
    printf(" %lld ns/ciclo\n", lat_henc);
    printf("  %-28s ", "Hash Open Addr. (otim.):");
    printBarra(lat_hotm, mx_lat, 28, CYN);
    printf(" %lld ns/ciclo\n", lat_hotm);

    printLinha('=', 57);
    printf(DIM "  Complexidade teorica: Hash O(1) | AVL O(log n) | Lista O(n)\n" RST);

    // EXPORTA RESULTADOS PARA CSV (para graficos do relatorio)
    ofstream csv("resultados_benchmark.csv");
    if (csv.is_open())
    {
        csv << "operacao,estrutura,ns_por_op\n";
        csv << "busca,hash_encadeamento," << ns_henc << "\n";
        csv << "busca,hash_open_addr," << ns_hotm << "\n";
        csv << "busca,avl," << ns_avl << "\n";
        csv << "busca,lista," << ns_lst << "\n";
        csv << "busca,bloom," << ns_blm << "\n";
        csv << "insercao,lista," << ins_lst << "\n";
        csv << "insercao,avl," << ins_avl << "\n";
        csv << "insercao,hash_encadeamento," << ins_hc << "\n";
        csv << "insercao,hash_open_addr," << ins_ho << "\n";
        csv << "insercao,kd_tree," << ins_kd << "\n";
        csv << "insercao,bloom," << ins_bl << "\n";
        csv << "remocao,lista," << rem_lst << "\n";
        csv << "remocao,avl," << rem_avl << "\n";
        csv << "remocao,hash_encadeamento," << rem_hc << "\n";
        csv << "remocao,hash_open_addr," << rem_ho << "\n";
        csv << "latencia_combinada,hash_encadeamento," << lat_henc << "\n";
        csv << "latencia_combinada,hash_open_addr," << lat_hotm << "\n";
        csv << "memoria_kb,hash_encadeamento," << mem_henc / 1024 << "\n";
        csv << "memoria_kb,hash_open_addr," << mem_hotm / 1024 << "\n";
        csv << "memoria_kb,bloom," << mem_blm / 1024 << "\n";
        csv << "memoria_kb,avl," << mem_avl / 1024 << "\n";
        if (rss >= 0)
            csv << "memoria_kb,rss_processo_real," << rss << "\n";
        csv.close();
        printf(GRN "  [CSV] Resultados salvos em 'resultados_benchmark.csv'\n" RST);
    }
}

// =============================================================================
//  12. TESTE DE ESCALABILIDADE (O(n) vs O(log n) vs O(1))
// =============================================================================

void testeEscalabilidade()
{
    printf("\n" CYN BOLD "  ESCALABILIDADE - Tempo medio de Busca por N\n" RST);
    printLinha('=', 57);
    printf(BOLD "  %-8s  %14s %14s %14s\n" RST,
           "N", "Lista O(n) ns", "AVL O(lgN) ns", "Hash O(1) ns");
    printLinha('-', 57);

    int tamanhos[] = {1000, 5000, 10000, 50000, 100000};
    const int NTAM = 5;
    double ns_lst[NTAM], ns_avl[NTAM], ns_hsh[NTAM];

    for (int t = 0; t < NTAM; t++)
    {
        int N = tamanhos[t];
        NoLista *lst = NULL;
        NoAVL *avl_t = NULL;
        TabelaHash *hsh = criarTabelaHash();
        EstadoAVLGuard guard;

        mt19937 gen(t * 31 + 7);
        uniform_real_distribution<double> d(0.001, 0.999);

        for (int i = 1; i <= N; i++)
        {
            double a = d(gen);
            lst = inserirLista(lst, a, i);
            avl_t = inserirAVL(avl_t, a, i, "esc");
            inserirHash(hsh, i, a);
        }

        // Coleta 32 pares (area, id) distribuidos para variar a busca entre iteracoes.
        const int NP = 32;
        double pa[NP];
        unsigned int pi_arr[NP];
        {
            int stp = max(1, N / NP);
            for (int p = 0; p < NP; p++)
            {
                pi_arr[p] = (unsigned int)(1 + p * stp);
                pa[p] = 0.0;
                buscarHash(hsh, pi_arr[p], &pa[p]);
            }
        }

        // Repeticoes ADAPTATIVAS: a lista e O(n) (cara), AVL/Hash sao baratas.
        // A lista no N pequeno e rapida: precisa de muitas repeticoes p/ superar
        // a resolucao do timer (~0.6ms); no N grande, poucas bastam.
        const int REP_LENTO = max(300, 10000000 / N);
        const int REP_RAPIDO = 2000000;
        // volatile counter impede o compilador de hoistar chamadas puras
        auto medirBloco = [&](int reps, auto fn) -> double
        {
            volatile int ctr = 0;
            auto t1 = chrono::high_resolution_clock::now();
            for (int i = 0; i < reps; i++)
            {
                int idx = (ctr++) % NP;
                fn(pa[idx], pi_arr[idx]);
            }
            long long tot = chrono::duration_cast<chrono::nanoseconds>(
                                chrono::high_resolution_clock::now() - t1)
                                .count();
            return (double)tot / reps;
        };

        ns_lst[t] = medirBloco(REP_LENTO, [&](double a, unsigned int)
                               { g_sink_b = buscarLista(lst, a); });
        ns_avl[t] = medirBloco(REP_RAPIDO, [&](double a, unsigned int)
                               { g_sink_b = (buscarAVL(avl_t, a) != NULL); });
        double ar = 0;
        ns_hsh[t] = medirBloco(REP_RAPIDO, [&](double, unsigned int id)
                               { g_sink_b = buscarHash(hsh, id, &ar); g_sink_d = ar; });

        // exibe ns/op com 2 decimais (revela valores sub-nanosegundo sem truncar p/ 0)
        printf("  %-8d  %11.2f ns %11.2f ns %11.2f ns\n",
               N, ns_lst[t], ns_avl[t], ns_hsh[t]);

        liberarLista(lst);
        liberarAVL(avl_t);
        liberarTabelaHash(hsh);
        // guard restaura g_total_avl_nodes automaticamente aqui
    }

    // Grafico proporcional para N=100k
    double mx = max({ns_lst[4], ns_avl[4], ns_hsh[4], 0.01});

    printf(BOLD "\n  Comparacao N=100.000 (proporcional ao mais lento):\n" RST);
    printf("  Lista  ");
    printBarra(ns_lst[4], mx, 34, RED);
    printf(" %.2f ns\n", ns_lst[4]);
    printf("  AVL    ");
    printBarra(max(ns_avl[4], 0.01), mx, 34, YEL);
    printf(" %.2f ns\n", ns_avl[4]);
    printf("  Hash   ");
    printBarra(max(ns_hsh[4], 0.01), mx, 34, GRN);
    printf(" %.2f ns\n", ns_hsh[4]);

    // Razoes de crescimento (usa N=1k e N=100k)
    printf(BOLD "\n  Razao de crescimento N=1k -> N=100k:\n" RST);
    double r_lst = ns_lst[0] > 0 ? ns_lst[4] / ns_lst[0] : 0;
    double r_avl = ns_avl[0] > 0 ? ns_avl[4] / ns_avl[0] : 0;
    double r_hsh = ns_hsh[0] > 0 ? ns_hsh[4] / ns_hsh[0] : 0;
    printf("  Lista: " RED "%.1fx" RST " (esperado ~100x para O(n))\n", r_lst);
    printf("  AVL:   " YEL "%.1fx" RST " (esperado ~%.1fx para O(log n))\n",
           r_avl, log2(100000.0) / log2(1000.0));
    printf("  Hash:  " GRN "%.1fx" RST " (esperado ~1x  para O(1))\n", r_hsh);
    printf(DIM "  A Lista cresce ALEM de 100x por efeito de cache: listas pequenas\n"
               "  cabem em cache (rapido por no); listas grandes sofrem cache miss\n"
               "  (lento por no). AVL/Hash ficam sub-microsegundo (cache L1).\n" RST);

    // EXPORTA para CSV (linha por N -> ideal para grafico de linha no relatorio)
    ofstream csv("resultados_escalabilidade.csv");
    if (csv.is_open())
    {
        csv << "N,lista_ns,avl_ns,hash_ns\n";
        for (int t = 0; t < NTAM; t++)
            csv << tamanhos[t] << "," << ns_lst[t] << "," << ns_avl[t] << "," << ns_hsh[t] << "\n";
        csv.close();
        printf(GRN "  [CSV] Dados salvos em 'resultados_escalabilidade.csv'\n" RST);
    }

    printLinha('=', 57);
    printf(DIM "  Busca no meio da lista = pior caso para O(n), melhor caso para O(log n)/O(1)\n" RST);
}

// =============================================================================
//  12b. DEMONSTRACAO DA RESTRICAO R2 (teto de memoria da AVL)
// =============================================================================

void testarRestricaoR2()
{
    printf("\n" CYN BOLD "  RESTRICAO R2 - Teto de Memoria da AVL\n" RST);
    printLinha('=', 57);
    printf("  O PDF (categoria 1) exige limitar o tamanho de uma estrutura.\n");
    printf("  Aqui a AVL tem um teto rigido: insercoes alem do limite sao\n");
    printf("  rejeitadas para nao ultrapassar o orcamento de memoria.\n");
    printLinha('-', 57);

    // Demonstra em uma AVL ISOLADA (nao mexe na arvore principal)
    const int TETO_DEMO = 5000;
    const int TENTATIVAS = 8000;

    // O guard salva o estado real da AVL e o restaura ao fim da funcao (RAII)
    EstadoAVLGuard guard;
    int teto_real = guard.limite;
    g_total_avl_nodes = 0;
    g_limite_avl = TETO_DEMO;
    g_avl_rejeitados = 0;

    printf("  Teto configurado:     " YEL "%d nos\n" RST, TETO_DEMO);
    printf("  Insercoes tentadas:   " YEL "%d\n" RST, TENTATIVAS);

    NoAVL *demo = NULL;
    mt19937 gen(123);
    uniform_real_distribution<double> d(0.001, 0.999);
    for (int i = 1; i <= TENTATIVAS; i++)
        demo = inserirAVL(demo, d(gen), i, "r2_demo");

    int aceitos = g_total_avl_nodes;
    long long rej = g_avl_rejeitados;

    printLinha('-', 57);
    printf("  Nos aceitos:          " GRN "%d" RST "  ", aceitos);
    printBarra(aceitos, TENTATIVAS, 24, GRN);
    printf("\n");
    printf("  Insercoes rejeitadas: " RED "%lld" RST "  ", rej);
    printBarra((double)rej, TENTATIVAS, 24, RED);
    printf("\n");
    printLinha('-', 57);

    if (aceitos == TETO_DEMO && rej == TENTATIVAS - TETO_DEMO)
        printf(GRN "  [OK] Restricao R2 funcionando: a AVL travou exatamente em %d nos.\n" RST, TETO_DEMO);
    else
        printf(YEL "  [!] Resultado inesperado (verifique a logica do teto).\n" RST);

    liberarAVL(demo);

    printLinha('=', 57);
    printf(DIM "  A arvore principal nao foi afetada (teste isolado). Teto real: %d nos.\n" RST,
           teto_real);
}


// ----- (origem: main.cpp) -----

#include <cstdio>
#include <iostream>
using namespace std;

// Definicao dos globais declarados em comum.h
int g_total_avl_nodes = 0;
bool g_janela_aberta = true;
int g_limite_avl = LIMITE_MAX_AVL;
long long g_avl_rejeitados = 0;
volatile bool g_sink_b = false;
volatile double g_sink_d = 0.0;

// =============================================================================
//  13. MENU PRINCIPAL
// =============================================================================

int main()
{
    habilitarANSI();

    NoLista *lista = NULL;
    TabelaHashOtm *hotm = criarHashOtm();
    NoAVL *avl = NULL;
    TabelaHash *hash = criarTabelaHash();
    NoKD *kd = NULL;
    BloomFilter *bloom = criarBloom();

    int op = 0;
    do
    {
        printf("\n");
        printHeader();
        printf(GRN "  >> DEMONSTRACAO: execute as opcoes 1 a 8 NESTA ORDEM <<\n" RST);
        printf(DIM "  --- 1) Carregar e tratar os dados -------------------\n" RST);
        printf(WHT "   1." RST " Ingerir Dataset            " DIM "(R11 rede, R18 sensor)\n" RST);
        printf(WHT "   2." RST " Histograma de Distribuicao " DIM "(ver dados brutos)\n" RST);
        printf(WHT "   3." RST " Filtrar Outliers Z/IQR     " DIM "(R2 limite, R18 poda)\n" RST);
        printf(DIM "  --- 2) Analisar -------------------------------------\n" RST);
        printf(WHT "   4." RST " Analise Espacial KD-Tree   " DIM "(produtividade por frame)\n" RST);
        printf(DIM "  --- 3) Medir desempenho -----------------------------\n" RST);
        printf(WHT "   5." RST " Teste de Escalabilidade    " DIM "(+CSV p/ relatorio)\n" RST);
        printf(WHT "   6." RST " Benchmark Completo         " DIM "(+CSV p/ relatorio)\n" RST);
        printf(DIM "  --- 4) Demonstrar restricoes ------------------------\n" RST);
        printf(WHT "   7." RST " Teste do Teto da AVL       " DIM "(restricao R2 memoria)\n" RST);
        printf(WHT "   8." RST " Teste de Bloqueio de CPU   " DIM "(restricao R9)\n" RST);
        printf(DIM "  --- Operacoes avulsas (opcional) --------------------\n" RST);
        printf(WHT "   9." RST " Buscar registro por ID\n");
        printf(WHT "  10." RST " Inserir registro individual\n");
        printf(WHT "  11." RST " Remover registro por ID\n");
        printf(DIM "  -----------------------------------------------------\n" RST);
        printf(WHT "   0." RST " Encerrar e liberar memoria\n");
        printLinha('-');
        printf(CYN "  Escolha: " RST);

        if (!(cin >> op))
        {
            cin.clear();
            cin.ignore(10000, '\n');
            op = -1;
        }

        switch (op)
        {

        case 1:
            carregarDataset(&lista, hotm, &avl, hash, &kd, bloom);
            printf(YEL "  Status AVL (R2): " RST "%d / %d nos (%.1f%% do teto)",
                   g_total_avl_nodes, g_limite_avl,
                   100.0 * g_total_avl_nodes / g_limite_avl);
            if (g_avl_rejeitados > 0)
                printf(RED "  | %lld rejeitadas pelo teto!\n" RST, g_avl_rejeitados);
            else
                printf(DIM "  | teto nao atingido (use a opcao 7 p/ ver R2 ativo)\n" RST);
            break;

        case 4:
            if (!kd)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }
            analiseEspacialKD(kd);
            break;

        case 3:
        {
            if (!lista)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }

            vector<double> vals;
            for (NoLista *p = lista; p; p = p->prox)
                vals.push_back(p->bounding_area);
            ResumoEstatistico r = calcularResumo(vals);

            // Conta quantas anomalias cada metodo capturaria
            auto contarFora = [&](double lo, double hi)
            {
                int c = 0;
                for (double v : vals)
                    if (v < lo || v > hi)
                        c++;
                return c;
            };
            int n_zsc = contarFora(r.zsc_inf, r.zsc_sup);
            int n_iqr = contarFora(r.iqr_inf, r.iqr_sup);

            printf("\n" CYN BOLD "  DETECCAO DE OUTLIERS - Z-SCORE vs IQR\n" RST);
            printLinha('-', 57);
            printf("  Z-Score (media +/- 3 sigma):\n");
            printf("    Faixa: " YEL "[%.4f , %.4f]" RST "  capturaria " RED "%d" RST " outliers\n",
                   r.zsc_inf, r.zsc_sup, n_zsc);
            printf("    " DIM "Problema: as anomalias inflam sigma (%.4f), alargando a faixa.\n" RST, r.desvio);
            printf("  IQR / Tukey (Q3 + 1.5*IQR):\n");
            printf("    Faixa: " YEL "[%.4f , %.4f]" RST "  capturaria " GRN "%d" RST " outliers\n",
                   r.iqr_inf, r.iqr_sup, n_iqr);
            printf("    " DIM "Robusto: baseado em mediana/quartis, imune a contaminacao.\n" RST);
            printLinha('-', 57);

            // Aplica o metodo robusto (IQR) na AVL
            printf("  Aplicando filtro IQR (robusto) na AVL (R18)...\n");
            int antes = g_total_avl_nodes;
            avl = podarAVL(avl, r.iqr_inf, r.iqr_sup);
            printf(GRN "  [OK] %d outliers removidos. Nos validos na AVL: %d\n" RST,
                   antes - g_total_avl_nodes, g_total_avl_nodes);
            printf(DIM "  Compare: o IQR removeu mais anomalias R18 que o z-score teria removido.\n" RST);
            printLinha('-', 57);
            break;
        }

        case 6:
            if (!lista)
            {
                printf(RED "  [!] Ingira os dados primeiro (opcao 1).\n" RST);
                break;
            }
            executarBenchmark(lista, hotm, avl, hash, bloom, lista->sample_id);
            break;

        case 5:
            testeEscalabilidade();
            break;

        case 2:
            exibirHistograma(lista, avl);
            break;

        case 8:
        {
            double ar = 0.0;
            unsigned int id_t = lista ? lista->sample_id : 1050;
            printf("\n" CYN BOLD "  TESTE R9 -Disponibilidade da CPU\n" RST);
            printLinha('-', 55);
            g_janela_aberta = false;
            printf("  Passo 1: " RED "CPU 100%% ocupada" RST " -tentando buscar ID %u:\n", id_t);
            buscarHashRestrita(hash, id_t, &ar);
            g_janela_aberta = true;
            printf("\n  Passo 2: " GRN "CPU liberada" RST " -tentando buscar ID %u:\n", id_t);
            if (buscarHashRestrita(hash, id_t, &ar))
                printf(GRN "  [OK] ID %u encontrado. Area: %.6f\n" RST, id_t, ar);
            else
                printf(YEL "  [!] ID %u nao encontrado (sem dados carregados).\n" RST, id_t);
            printLinha('-', 55);
            break;
        }

        case 10:
        {
            unsigned int nid = (unsigned int)lerInteiro("  ID do novo registro: ");
            double narea = lerDouble("  Bounding area:       ");
            lista = inserirLista(lista, narea, nid);
            inserirHashOtm(hotm, nid, narea);
            avl = inserirAVL(avl, narea, nid, "manual.png");
            inserirHash(hash, nid, narea);
            inserirBloom(bloom, nid);
            printf(GRN "  [OK] ID %u inserido em todas as estruturas.\n" RST, nid);
            break;
        }

        case 9:
        {
            unsigned int bid = (unsigned int)lerInteiro("  ID para buscar: ");
            double bar = 0.0;
            bool bBlm = buscarBloom(bloom, bid);
            bool bHsh = buscarHash(hash, bid, &bar);
            double bOptArea = 0.0;
            bool bOpt = buscarHashOtm(hotm, bid, &bOptArea);
            NoAVL *bAvl = buscarAVL(avl, bar);

            printf("\n" CYN BOLD "  RESULTADO -ID %u\n" RST, bid);
            printLinha('-', 55);
            // Bloom nao suporta delecao: pode retornar true mesmo apos remocao (falso positivo)
            if (bBlm && !bHsh)
                printf("  Bloom Filter:    " YEL "FALSO POSITIVO (elemento foi removido, Bloom nao suporta delecao)\n" RST);
            else if (bBlm)
                printf("  Bloom Filter:    " GRN "Possivel membro (confirmar nas estruturas exatas abaixo)\n" RST);
            else
                printf("  Bloom Filter:    " RED "AUSENTE com certeza (nenhum falso negativo possivel)\n" RST);
            printf("  Hash Encadeada:  %s  area=%.6f\n", bHsh ? GRN "ENCONTRADO" RST : YEL "nao encontrado" RST, bar);
            printf("  Hash Otimizada:  %s  area=%.6f\n", bOpt ? GRN "ENCONTRADO" RST : YEL "nao encontrado" RST, bOptArea);
            printf("  AVL:             %s\n", bAvl ? GRN "ENCONTRADO" RST : YEL "nao encontrado (pode ter sido podado)" RST);
            printLinha('-', 55);
            break;
        }

        case 11:
        {
            unsigned int rid = (unsigned int)lerInteiro("  ID para remover: ");
            double rar = 0.0;
            if (!buscarHash(hash, rid, &rar))
            {
                printf(YEL "  [!] ID %u nao encontrado.\n" RST, rid);
                break;
            }
            lista = removerLista(lista, rid);
            removerHashOtm(hotm, rid);
            avl = removerAVL(avl, rar, rid);
            removerHash(hash, rid);
            printf(GRN "  [OK] ID %u removido de todas as estruturas.\n" RST, rid);
            break;
        }

        case 7:
            testarRestricaoR2();
            break;

        case 0:
            printf(CYN "  Desalocando memoria e encerrando...\n" RST);
            break;

        default:
            printf(RED "  [!] Opcao invalida.\n" RST);
        }
    } while (op != 0);

    liberarLista(lista);
    liberarHashOtm(hotm);
    liberarAVL(avl);
    liberarTabelaHash(hash);
    liberarKD(kd);
    liberarBloom(bloom);
    return 0;
}
