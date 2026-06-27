#include "estruturas.h"
#include "comum.h"
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

